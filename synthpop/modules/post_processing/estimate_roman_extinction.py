"""
Postprocessing module to improve the estimation of F146 wide-filter extinction
with a polynomial function of A_Ks and narrower filter colors.
"""

__all__ = ["EstimateRomanExtinction", ]
__author__ = "M.J. Huston"
__date__ = "2026-04-26"

import pandas as pd 
import numpy as np
from ._post_processing import PostProcessing
from .convert_mist_mags import ConvertMistMags
import os.path
import json
import warnings
import pdb
from ...synthpop_utils.utils_functions import combine_system_mags, get_primary_mags


class EstimateRomanExtinction(PostProcessing):
    """
    Postprocessing module to apply a correction to F146 extinction
    
    Attributes
    ----------
    """

    def __init__(self, model, logger, mag_sys='None', use_low_extinction=False, **kwargs):
        super().__init__(model,logger, **kwargs)
        if mag_sys != 'AB':
            raise ValueError('To use CorrectF146Extinction, magnitudes must be in AB. Confirm that you '
                'have converted to AB mags before using this module, then set this module\'s mag_sys ' 
                'kwarg to \'AB\' before running. You can put the ConvertMistMags module before this '
                'one in your post_processing_kwargs to do the conversion.')

        # Load up the Roman fit results
        self.use_low_extinction = use_low_extinction
        current_dir = os.path.dirname(os.path.abspath(__file__))
        with open(f'{current_dir}/roman_fits_abs_AKs5.json', 'r') as f:
            self.roman_ext_fits = json.load(f)
        with open(f'{current_dir}/roman_fits_abs_AKs1.json', 'r') as f:
            self.roman_ext_fits_lowext = json.load(f)
        self.roman_filter_list = ['f062', 'f087', 'f106', 'f129', 'f158', 'f184', 'f213', 'f146']

    @staticmethod
    def generic_extinction_polynomial(AKs_C, coeffs, order):
        """
        Generic polynomial function of flexible order. Any number of colors
        may be input, and cross-terms are computed between A_Ks and each color,
        not between different colors.
        
        Identical function to that used to run the fits.

        Parameters:
        -----------
        AKs_C : ndarray
            array of extinction and color values to compute extinction estimate
        coeffs : ndarray
            array of coefficients for the polynomial function
        order : int
            polynomial order for the extinction coefficient function

        Returns:
        --------
        ext_ests : ndarray
            extinction estimates for each star in the AKs_C table
        """
        n_colors = AKs_C.shape[1] - 1 
        n_terms_AKs = order + 1
        n_terms_per_color = order * (order + 1) // 2
        n_terms = n_terms_AKs + (n_colors * n_terms_per_color)
        assert n_terms == len(coeffs)

        AKs = AKs_C[:,0]
        var_terms = [AKs**p for p in range(order+1)]
        n_colors = AKs_C.shape[1]-1
        for i in range(0, n_colors):
            Ci = AKs_C[:, i+1]
            for p in range(1, order+1):
                for q in range(order+1):
                    if p + q <= order:
                        var_terms.append((Ci**p) * (AKs**q))
        terms_mat = np.column_stack(var_terms)
        val = terms_mat @ np.array(coeffs)
        return val * AKs

    def get_roman_extinction_sim(self, catalog):
        """
        Roman extinction estimator for simulations. Assumes all Roman filter photometry 
        is provided and in absolute AB mags.
        
        Parameters:
        -----------
        catalog : pd.DataFrame, astropy.table.Table, or similar
            required columns: A_Ks, f062, f087, f106, f129, f158, f184, f213, and f146
        low_extinction=False : boolean
            if all A_Ks<=1, use the alternate lower order correction. if any A_Ks>1, a warning
            will be printed, and the higher order correction will be used
            
        Returns:
        --------
        extinctions : dict
            entries of '<filter>':[<ext_star1>, <ext_star2>, ...] for each filter
        """
        # Select the appropriate fit_dict
        fit_dict = self.roman_ext_fits
        if self.use_low_extinction and np.all(catalog['A_Ks']<=1):
            fit_dict = self.roman_ext_fits_lowext
        elif self.use_low_extinction:
            warnings.warn("low_extinction set to True, but some A_Ks > 1. "
                          "switching to 0 <= A_Ks <= 5 fit.")
        
        # Iterate over the filters
        result = {}
        for filt in self.roman_filter_list:
            filt_valid = (filt in catalog)
            filt_fit = fit_dict[filt]
            colors = filt_fit['colors']
            coeffs = filt_fit['coefficients']
            order = filt_fit['order']
            print(f"estimating {filt} extinction using {colors} and order={order} function")

            columns = [catalog['A_Ks']]
            for c in colors:
                f1,f2 = c.split('_')
                filt_valid *= (f1 in catalog)
                filt_valid *= (f2 in catalog)
                if filt_valid:
                    columns.append(catalog[f1]-catalog[f2])
                else:
                    columns.append(np.ones(len(catalog))*np.nan)
            AKs_C = np.stack(columns,axis=1)
            
            ext_filt = self.generic_extinction_polynomial(AKs_C, coeffs, order)
            result['A_'+filt] = ext_filt
            
        return result

    def do_post_processing(self, systems: pd.DataFrame, companions: pd.DataFrame):
        """
        Run the process
        """
        if 'W146' in systems:
            if "2MASS_Ks" in systems:
                systems.loc[:,"K213"] = systems['2MASS_Ks'] + 1.834505
                companions.loc[:,"K213"] = companions['2MASS_Ks'] + 1.834505
                warnings.warn("K213 missing from MISTv1, estimating from 2MASS_Ks "
                    "assuming 2MASS mags are in Vegamag")
                self.model.parms.eff_wavelengths['K213'] = self.model.parms.eff_wavelengths['2MASS_Ks']
            mag_cols = ["R062", "Z087", "Y106", "J129", "H158", "F184", "K213", "W146"]
        elif 'm_roman_f146' in systems:
            mag_cols = [f'm_roman_{f}' for f in self.roman_filter_list]

        # Get and convert extinction to A_Ks if needed
        A_Ks = systems[self.model.populations[0].extinction.A_or_E_type].to_numpy()
        if self.model.populations[0].extinction.A_or_E_type != 'A_Ks':
            A_Ks *= self.model.populations[0].extinction.Alambda_Amap(2.152152)
            systems.loc[:,'A_Ks'] = A_Ks
        ext_est_dict = {"A_Ks":A_Ks}

        # Do systems in normal case
        if (companions is None) or (not self.model.parms.combine_system_mags):
            # Get absolute mags if needed
            if self.model.parms.obsmag:
                # Convert needed mags back to absolute
                for i,f in enumerate(mag_cols):
                    if f in systems:
                        ext_est_dict[self.roman_filter_list[i]] = (systems[f].to_numpy() - 
                                5*np.log10(100*systems['Dist'].to_numpy()) - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[f]) * \
                                systems[self.model.populations[0].extinction.A_or_E_type].to_numpy())
                    else:
                        warnings.warn(f"Column {f} not found in systems, some extinction corrections "
                            "will not be estimated.")
            else:
                for i,f in enumerate(mag_cols):
                    if f in systems:
                        ext_est_dict[self.roman_filter_list[i]] = systems[f].to_numpy()
                    else:
                        warnings.warn(f"Column {f} not found in systems, some extinction corrections "
                            "will not be estimated.")

            # Calculate extinction in the filters from the A_Ks and absolute colors
            ext_ests = self.get_roman_extinction_sim(ext_est_dict)
            for i,f in enumerate(self.roman_filter_list):
                sp_f = mag_cols[i]
                if not np.all(np.isnan(ext_ests['A_'+f])):
                    systems.loc[:,f"A_{sp_f}"] = ext_ests['A_'+f]
                    if self.model.parms.obsmag:
                        systems.loc[:,sp_f] += ext_ests['A_'+f] - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[sp_f]) * \
                                systems[self.model.populations[0].extinction.A_or_E_type].to_numpy()
                else:
                    warnings.warn(f"{sp_f} extinction correction could not be estimated.")

        # Do companions
        if companions is not None:
            comps_abs_mags = companions[['system_idx']].copy()
            A_Ks_series = pd.Series(A_Ks, index=systems['system_idx'].to_numpy())
            ext_est_dict['A_Ks'] = A_Ks_series[companions['system_idx'].to_numpy()]
            Dist_series = pd.Series(systems['Dist'].to_numpy(), index=systems['system_idx'].to_numpy())
            # Get absolute mags if needed
            if self.model.parms.obsmag:
                # Convert needed mags back to absolute
                for i,f in enumerate(mag_cols):
                    if f in companions:
                        abs_mag_f = (companions[f].to_numpy() - 
                                5*np.log10(100*Dist_series[companions['system_idx'].to_numpy()].to_numpy()) - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[f]) * \
                                A_Ks_series[companions['system_idx'].to_numpy()].to_numpy())
                        ext_est_dict[self.roman_filter_list[i]] = abs_mag_f
                        comps_abs_mags.loc[:,self.roman_filter_list[i]] = abs_mag_f
            else:
                for i,f in enumerate(mag_cols):
                    if f in companions:
                        abs_mag_f = companions[f].to_numpy()
                        ext_est_dict[self.roman_filter_list[i]] = abs_mag_f
                        comps_abs_mags.loc[:,self.roman_filter_list[i]] = abs_mag_f

            ext_ests = self.get_roman_extinction_sim(ext_est_dict)
            for i,f in enumerate(self.roman_filter_list):
                sp_f = mag_cols[i]
                if not np.all(np.isnan(ext_ests['A_'+f])):
                    companions.loc[:,f"A_{sp_f}"] = ext_ests['A_'+f]
                    if self.model.parms.obsmag:
                        companions.loc[:,sp_f] += ext_ests['A_'+f] - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[sp_f]) * \
                                A_Ks_series[companions['system_idx'].to_numpy()].to_numpy()

        # Do combined mag case systems table
        if (companions is not None) and (self.model.parms.combine_system_mags):
            combo_abs_mags = systems[['system_idx','n_companions', 'A_Ks']].copy()
            if self.model.parms.combine_system_mags:
                prim_abs_mags = systems[['system_idx','n_companions', 'A_Ks']].copy()
            # Get absolute mags if needed
            if self.model.parms.obsmag:
                # Convert needed mags back to absolute
                for i,f in enumerate(mag_cols):
                    if f in systems:
                        abs_mag_f = (systems[f].to_numpy() - 
                                5*np.log10(100*systems['Dist'].to_numpy()) - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[f]) * \
                                systems[self.model.populations[0].extinction.A_or_E_type].to_numpy())
                        if self.model.parms.combine_system_mags:
                            combo_abs_mags.loc[:,self.roman_filter_list[i]] = abs_mag_f
                        else:
                            prim_abs_mags.loc[:,self.roman_filter_list[i]] = abs_mag_f
                    else:
                        warnings.warn(f"Column {f} not found in systems, some extinction corrections "
                            "will not be estimated.")
            else:
                for i,f in enumerate(mag_cols):
                    if f in systems:
                        abs_mag_f = systems[f].to_numpy()
                        if self.model.parms.combine_system_mags:
                            combo_abs_mags.loc[:,self.roman_filter_list[i]] = abs_mag_f
                        else:
                            prim_abs_mags.loc[:,self.roman_filter_list[i]] = abs_mag_f
                    else:
                        warnings.warn(f"Column {f} not found in systems, some extinction corrections "
                            "will not be estimated.")

            if self.model.parms.combine_system_mags:
                prim_abs_mags = get_primary_mags(combo_abs_mags, comps_abs_mags, self.roman_filter_list)
                prim_app_mags = systems[['system_idx']].copy()

            # Calculate extinction in the filters from the A_Ks and absolute colors
            ext_ests = self.get_roman_extinction_sim(prim_abs_mags)
            for i,f in enumerate(self.roman_filter_list):
                sp_f = mag_cols[i]
                if not np.all(np.isnan(ext_ests['A_'+f])):
                    systems.loc[:,f"A_{sp_f}"] = ext_ests['A_'+f]
                    # So, in the absolute mag case, we are good now.
                    # Next case: apparent mag case, systems not combined. 
                    # Here, just correct the primary star as you would a single star
                    if self.model.parms.obsmag and (not self.model.parms.combine_system_mags):
                        systems.loc[:,sp_f] += ext_ests['A_'+f] - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[sp_f]) * \
                                systems[self.model.populations[0].extinction.A_or_E_type].to_numpy()
                    # Next case: absolute mags for combined systems
                    elif self.model.parms.obsmag and self.model.parms.combine_system_mags:
                        prim_app_mags.loc[:,sp_f] = prim_abs_mags[self.roman_filter_list[i]].to_numpy() + \
                                                    ext_ests['A_'+f] + 5*np.log10(100*systems['Dist'].to_numpy())
                else:
                    warnings.warn(f"{sp_f} extinction correction could not be estimated.")
            if self.model.parms.obsmag and self.model.parms.combine_system_mags:
                system_app_mags = combine_system_mags(prim_app_mags, companions, mag_cols)
                for i,f in enumerate(self.roman_filter_list):
                    sp_f = mag_cols[i]
                    systems.loc[:,sp_f] = system_app_mags[sp_f].to_numpy()

        return systems, companions
