"""
Postprocessing module to improve the estimation extinction in photometric filters
with a polynomial function of A_Ks and absolute colors.
"""

__all__ = ["EstimateRomanExtinction", ]
__author__ = "M.J. Huston"
__date__ = "2026-04-26"

import pandas as pd 
import numpy as np
from ._post_processing import PostProcessing
import os.path
import json
import warnings
import pdb
from astropy.table import Table
from ...synthpop_utils.utils_functions import combine_system_mags, get_primary_mags
from ..evolution._evolution import EVOLUTION_DIR

class ExtinctionEstimator(PostProcessing):
    """
    Postprocessing module to apply a correction to F146 extinction
    
    Attributes
    ----------
    """

    def __init__(self, model, logger, **kwargs):
        super().__init__(model,logger, **kwargs)
        # Load up the fit results
        current_dir = os.path.dirname(os.path.abspath(__file__))
        ext_cor_tab = Table.read(f'{current_dir}/extinction_correction_table.dat', format='ascii.mrt')
        ext_cor_tab = ext_cor_tab.filled(999)
        self.fit_dict = {}
        color_cols = [col for col in ext_cor_tab.colnames if col.startswith('Color')]
        a_cols = [col for col in ext_cor_tab.colnames if col.startswith('a_')]
        coeff_cols = [col for col in ext_cor_tab.colnames if (col[:2] in ['a_', 'b_'])]
        self.fit_order = len(a_cols)-1
        for i in range(len(ext_cor_tab)):
            filt = ext_cor_tab['Filter'][i]
            self.fit_dict[filt] = {}
            self.fit_dict[filt]['colors'] = [ext_cor_tab[col][i] for col in color_cols if (ext_cor_tab[col][i]!="999")]
            self.fit_dict[filt]['order'] = self.fit_order
            self.fit_dict[filt]['coefficients'] = [ext_cor_tab[col][i] for col in coeff_cols if ~(ext_cor_tab[col][i]==999.0)]
        with open(f"{EVOLUTION_DIR}/spisea_photometric_system_conversions.json") as f:
            self.photsys_convert = json.load(f)

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
        # fit_dict = self.roman_ext_fits
        # if self.use_low_extinction and np.all(catalog['A_Ks']<=1):
        #     fit_dict = self.roman_ext_fits_lowext
        # elif self.use_low_extinction:
        #     warnings.warn("low_extinction set to True, but some A_Ks > 1. "
        #                   "switching to 0 <= A_Ks <= 5 fit.")
        
        # Iterate over the filters
        catalog = self.convert_mags_to_ab(catalog)
        result = {}
        for filt in self.filter_list:
            filt_valid = (filt in catalog)
            filt_fit = self.fit_dict[filt]
            colors = filt_fit['colors']
            coeffs = filt_fit['coefficients']
            order = filt_fit['order']
            self.logger.info(f"Estimating {filt} extinction using {colors} and order={order} function")

            columns = [catalog['A_Ks']]
            for c in colors:
                f1,f2,_ = c.split('-')
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

    def convert_mags_to_ab(self, ext_est_dict):
        for i, f in enumerate(self.full_filter_list_obs_str):
            if self.model.parms.photsys_dict[self.model.parms.bands[i]] == 'AB':
                pass 
            elif self.model.parms.photsys_dict[self.model.parms.bands[i]] == 'Vega':
                ext_est_dict[self.full_filter_list[i]] += self.photsys_convert['AB'][f]
            elif self.model.parms.photsys_dict[self.model.parms.bands[i]] == 'ST':
                ext_est_dict[self.full_filter_list[i]] -= self.photsys_convert['ST'][f]
                ext_est_dict[self.full_filter_list[i]] += self.photsys_convert['AB'][f]
        return ext_est_dict 

    def do_post_processing(self, systems: pd.DataFrame, companions: pd.DataFrame):
        """
        Run the process
        """

        # Catch case where we don't have K213 and need to swap in 2MASS_Ks
        if 'W146' in systems:
            if "2MASS_Ks" in systems:
                systems.loc[:,"K213"] = systems['2MASS_Ks']
                companions.loc[:,"K213"] = companions['2MASS_Ks']
                warnings.warn("K213 missing from MISTv1, estimating from 2MASS_Ks.")
                self.model.parms.eff_wavelengths['K213'] = self.model.parms.eff_wavelengths['2MASS_Ks']
                self.model.parms.photsys_dict["K213"] = self.model.parms.photsys_dict["2MASS_Ks"]
                self.model.parms.bands += ["K213"]

        # Set up filter sets
        if self.model.parms.star_generator=='SpiseaGenerator':
            from spisea.synthetic import get_obs_str
            self.full_filter_list_obs_str = [get_obs_str(f) for f in self.model.parms.bands]
        else:
            from ..evolution.mist import get_spisea_obs_str
            self.full_filter_list_obs_str = [get_spisea_obs_str(f)  for f in self.model.parms.bands]
        self.full_filter_list = [f.replace(',','_')  for f in self.full_filter_list_obs_str]
        self.correct_mag_cols = self.model.parms.bands
        self.filter_list = [f for f in self.full_filter_list if (f in self.fit_dict.keys())]
        if len(self.filter_list)<len(self.full_filter_list):
            not_in_cor_list = [f for f in self.full_filter_list if (f not in self.filter_list)]
            warnings.warn(f"Filters {not_in_cor_list} do not have tabulated extinction"
                " corrections and will use the default effective wavelength extinction law scaling.")
            self.correct_mag_cols = [self.correct_mag_cols[i] for i in range(len(self.correct_mag_cols)) 
                                     if (self.full_filter_list[i] in self.filter_list)]

        # Get and convert extinction to A_Ks if needed
        A_Ks = systems[self.model.populations[0].extinction.A_or_E_type].to_numpy()
        if self.model.populations[0].extinction.A_or_E_type != 'A_Ks':
            A_Ks *= self.model.populations[0].extinction.Alambda_Amap(2.152152)
            systems.loc[:,'A_Ks'] = A_Ks
        ext_est_dict = {"A_Ks":A_Ks}

        # Do systems in normal case
        if (companions is None) or (not self.model.parms.combine_system_mags):
            # Get absolute mags for extinction estimator
            for i,f in enumerate(self.model.parms.bands):
                if self.model.parms.obsmag:
                    ext_est_dict[self.full_filter_list[i]] = (systems[f].to_numpy() - 
                        5*np.log10(100*systems['Dist'].to_numpy()) - \
                        self.model.populations[0].extinction.Alambda_Amap(
                            self.model.parms.eff_wavelengths[f]) * \
                        systems[self.model.populations[0].extinction.A_or_E_type].to_numpy())
                else:
                    ext_est_dict[self.full_filter_list[i]] = systems[f].to_numpy()

            # Calculate extinction in the filters from the A_Ks and absolute colors
            ext_ests = self.get_roman_extinction_sim(ext_est_dict)
            for i,f in enumerate(self.filter_list):
                sp_f = self.correct_mag_cols[i]
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
            # Convert needed mags back to absolute
            for i,f in enumerate(self.model.parms.bands):
                if f in companions:
                    if self.model.parms.obsmag:
                        abs_mag_f = (companions[f].to_numpy() - 
                                5*np.log10(100*Dist_series[companions['system_idx'].to_numpy()].to_numpy()) - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[f]) * \
                                A_Ks_series[companions['system_idx'].to_numpy()].to_numpy())
                    else:
                        abs_mag_f = companions[f].to_numpy()
                    ext_est_dict[self.full_filter_list[i]] = abs_mag_f
                    comps_abs_mags.loc[:,self.full_filter_list[i]] = abs_mag_f

            ext_ests = self.get_roman_extinction_sim(ext_est_dict)
            for i,f in enumerate(self.filter_list):
                sp_f = self.correct_mag_cols[i]
                if not np.all(np.isnan(ext_ests['A_'+f])):
                    companions.loc[:,f"A_{sp_f}"] = ext_ests['A_'+f]
                    if self.model.parms.obsmag:
                        companions.loc[:,sp_f] += ext_ests['A_'+f] - \
                                self.model.populations[0].extinction.Alambda_Amap(
                                    self.model.parms.eff_wavelengths[sp_f]) * \
                                A_Ks_series[companions['system_idx'].to_numpy()].to_numpy()

        # Do combined mag case systems table
        if (companions is not None) and (self.model.parms.combine_system_mags):
            prim_abs_mags = systems[['system_idx','n_companions', 'A_Ks']].copy()
            # Get absolute mags, calculating if needed
            for i,f in enumerate(self.model.parms.bands):
                if self.model.parms.obsmag:
                    abs_mag_f = (systems[f].to_numpy() - 
                        5*np.log10(100*systems['Dist'].to_numpy()) - \
                        self.model.populations[0].extinction.Alambda_Amap(
                            self.model.parms.eff_wavelengths[f]) * \
                        systems[self.model.populations[0].extinction.A_or_E_type].to_numpy())
                else:
                    abs_mag_f = systems[f].to_numpy()
                prim_abs_mags.loc[:,self.full_filter_list[i]] = abs_mag_f
            # Subtract the companion mags if needed
            if self.model.parms.combine_system_mags:
                prim_abs_mags = get_primary_mags(prim_abs_mags, comps_abs_mags, self.full_filter_list)

            # Calculate extinction in the filters from the A_Ks and absolute colors
            ext_ests = self.get_roman_extinction_sim(prim_abs_mags)

            # Put the values in the table, and re-calculate apparent mags if relevant
            for i,f in enumerate(self.filter_list):
                sp_f = self.correct_mag_cols[i]
                if not np.all(np.isnan(ext_ests['A_'+f])):
                    systems.loc[:,f"A_{sp_f}"] = ext_ests['A_'+f]
                    # So, in the absolute mag case, we are good now.
                    # Apparent mag case: calculate primary apparent mags first, add companions later if needed
                    if self.model.parms.obsmag:
                        systems.loc[:,sp_f] = prim_abs_mags[self.filter_list[i]].to_numpy() + \
                                                    ext_ests['A_'+f] + 5*np.log10(100*systems['Dist'].to_numpy())
                else:
                    warnings.warn(f"{sp_f} extinction correction could not be estimated.")
            if self.model.parms.obsmag and self.model.parms.combine_system_mags:
                systems = combine_system_mags(systems, companions, self.correct_mag_cols)

        return systems, companions
