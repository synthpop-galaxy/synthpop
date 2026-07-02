"""
Post-processing to convert the output into the PopSyCLE input format,
ready to be plugged in at the calc_events stage.
Note: obsmag must be set to FALSE in config file to use this module
"""

__author__ = "M.J. Huston, S. Brooke, A. Kim"

from ._post_processing import PostProcessing
import time
import pandas as pd
import numpy as np
import h5py
import os
from popsycle.synthetic import _get_bin_edges, _bin_lb_hdf5
import pdb

filter_set_dict = {'ubv': ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']}

filter_matching_mist = {'ubv_J': "UKIDSS_J",
                        'ubv_H': "UKIDSS_H",
                        'ubv_K': "UKIDSS_K",
                        'ubv_U': "Bessell_U",
                        'ubv_I': "Bessell_I",
                        'ubv_B': "Bessell_B",
                        'ubv_V': "Bessell_V",
                        'ubv_R': "Bessell_R"
                       }

filter_spisea_dict = {'UKIDSS_J': "m_ukirt_J",
                        'UKIDSS_H': "m_ukirt_H",
                        'UKIDSS_K': "m_ukirt_K",
                        'Bessell_U': "m_ubv_U",
                        'Bessell_I': "m_ubv_I",
                        'Bessell_B': "m_ubv_B",
                        'Bessell_V': "m_ubv_V",
                        'Bessell_R': "m_ubv_R",
                     }

roman_filters_dict = {'roman_f062': "m_roman_f062",
                      'roman_f087': "m_roman_f087",
                      'roman_f106': "m_roman_f106",
                      'roman_f129': "m_roman_f129",
                      'roman_f158': "m_roman_f158",
                      'roman_f146': "m_roman_f146",
                      'roman_f184': "m_roman_f184",
                      'roman_f213': "m_roman_f213"
                     }

synthpop_nonmag_cols = ['l', 'b', 'Dist',
                        'x', 'y', 'z',
                        'vr_bc', 'mul','mub',
                        'U', 'V', 'W',
                        'iMass','Mass',
                        'log_L', 'log_g', 'log_Teff', 'Fe/H_initial', 'age',
                        'pop', 'phase', 'n_companions', 'system_Mass'
                        ]

synthpop_nonmag_bin_cols = ['system_idx', 'period', 'eccentricity', '2MASS_Ks',
                           'star_mass', 'log_R', 'log_a']

popsycle_nonmag_bin_cols = ['system_idx', 'zams_mass', 'Teff', 'L',
                            'logg', 'isWR', 'mass', 'phase', 'metallicity',
                            'log_a', 'e', 'i', 'Omega', 'omega', 'zams_mass_prim',
                            'glat','glon'
                            ]

popsycle_nonmag_cols = ['glat', 'glon', 'rad',
                        'px', 'py', 'pz', 
                        'vr', 'mu_lcosb', 'mu_b', 
                        'vx', 'vy', 'vz', 
                        'zams_mass', 'mass', 'systemMass', 
                        'mbol', 'grav', 'Teff', 'L', 'feh', 'age',
                        'exbv', 'popid',
                        'isMultiple', 'N_companions', 'rem_id', 'obj_id'] 

class PopsyclePostProcessing(PostProcessing):

    def __init__(self, model, logger, bin_edges_number=None,
                 filter_sets=['ubv'], **kwargs):
        super().__init__(model, logger, **kwargs)
        self.bin_edges_number = bin_edges_number
        #self.filter_sets = filter_sets
        self.mag_cols = []
        self.synthpop_mag_cols = []
        for fset in filter_sets:
            self.mag_cols += [fset+'_'+f for f in filter_set_dict[fset]]
            # if model.parms.star_generator != 'SpiseaGenerator':
                # self.synthpop_mag_cols += [filter_spisea_dict[fset+'_'+f] for f in filter_set_dict[fset]]
            self.synthpop_mag_cols += [filter_matching_mist[fset+'_'+f] for f in filter_set_dict[fset]]
        self.star_generator = model.parms.star_generator
        if hasattr(model.parms, "binning_procedure"):
            self.binning_procedure = True
        else:
            self.binning_procedure = False
        self.multiplicity = model.parms.multiplicity_kwargs

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Converts DataFrame into format needed for input to PopSyCLE as a replacement
        for Galaxia, saving the file to the set file name + '_psc.h5'.
        """
        self.logger.info(f"Beginning PopSyCLE postprocessing.")

        if self.multiplicity == None:
            companion_df = pd.DataFrame(columns=synthpop_nonmag_bin_cols + synthpop_nonmag_cols)  # Empty Dataframe for operations if singles

        
        nan_sorted_columns = companion_df.isna().sum().sort_values(ascending=False).index  # Sorts column names by amount of NaNs
        companion_df = companion_df.reindex(columns=nan_sorted_columns)  # Resorts columns to put NaN filled columns at the end
        companion_df = companion_df.loc[:,~companion_df.columns.duplicated()].copy()  # Removes duplicated column names

        if self.star_generator == 'SpiseaGenerator':
            system_df.rename(columns={filter_spisea_dict[f]:f for f in filter_spisea_dict},
                         inplace=True, errors='ignore')  # Rename SPISEA filter columns to match MIST filters
            companion_df.rename(columns={filter_spisea_dict[f]:f for f in filter_spisea_dict},
                         inplace=True, errors='ignore')  # Rename SPISEA filter columns to match MIST filters
            system_df.rename(columns={roman_filters_dict[f]:f for f in roman_filters_dict},
                         inplace=True)

            
        # Drop unused data
        cols_to_cut = []
        for col in system_df.keys():
            if col not in (synthpop_nonmag_cols+self.synthpop_mag_cols+list(roman_filters_dict)+
                           self.mag_cols+synthpop_nonmag_bin_cols+
                           [self.model.populations[0].extinction.A_or_E_type]):
                cols_to_cut.append(col)

        system_df.drop(columns=cols_to_cut, inplace=True)
        self.logger.info("Unused columns dropped")
        self.output_root = f"{self.model.get_filename(self.model.l_deg, self.model.b_deg)}_psc"

        bin_cols_to_cut = [] # bandaid
        for col in companion_df.keys():
            if col not in (synthpop_nonmag_cols+self.synthpop_mag_cols+list(roman_filters_dict)+ self.mag_cols+synthpop_nonmag_bin_cols+ [self.model.populations[0].extinction.A_or_E_type]+list(filter_spisea_dict)):
                bin_cols_to_cut.append(col)
        companion_df.drop(columns=bin_cols_to_cut, inplace=True)
                               

        system_df['isMultiple'] = system_df['n_companions'] 
        system_df.loc[system_df['isMultiple'] != 0, 'isMultiple'] = 1  # isMultiple True where n_companions is not 0

        if system_df['isMultiple'].any():
            map = system_df.set_index('system_idx')['iMass'].squeeze()
            companion_df['zams_mass_prim'] = companion_df['system_idx'].map(map)  # Primary star mass tracked in companion dataframe

        # Translate extinction to Ebv extinction
        if not self.model.populations[0].extinction.A_or_E_type=="E(B-V)":
            extinction = self.model.populations[0].extinction
            extinction_type = self.model.populations[0].extinction.A_or_E_type
            system_df.rename(columns={extinction_type:'ext_orig'}, inplace=True)
            #ext_in_map = system_df[extinction_type].to_numpy()
            Av_Aref = extinction.extinction_at_lambda(0.544579, 1.0)
            Ab_Aref = extinction.extinction_at_lambda(0.438074, 1.0)
            pd.eval("exbv = (Ab_Aref - Av_Aref)*system_df.ext_orig", target=system_df, inplace=True)
            system_df.drop(columns=['ext_orig'], inplace=True)
        else:
            system_df.rename(columns={"E(B-V)":'exbv'}, inplace=True)
        
        # create log (with same info as galaxia log)
        #dtype = [('latitude', 'f8'), ('longitude', 'f8'), ('surveyArea', 'f8')]
        #log = np.zeros(1, dtype=dtype)
        latitude = self.model.l_deg
        longitude = self.model.b_deg
        surveyArea = self.model.field_scale
        if self.model.field_scale_unit=='sr':
            surveyArea *= (180/np.pi)**2
            
        os.makedirs('/'.join(self.output_root.split('/')[:-1]), exist_ok=True)
        if not self.binning_procedure:
            if os.path.exists(self.output_root + '_synthpop_params.txt'):
                with open(self.output_root + '_synthpop_params.txt', 'r') as params_file:
                    lines = params_file.read()
            else:
                lines = ""
            with open(self.output_root + '_synthpop_params.txt', 'w') as params_file:
                params_file.write(f"seed {self.model.parms.random_seed}\n")
                params_file.write(lines)
            self.logger.info("Parameter file written")

        system_df.rename(columns={'iMass': 'zams_mass','Mass': 'mass',
                                 'x': 'px', 'y': 'py', 'z': 'pz',
                                 'U': 'vx', 'V': 'vy', 'W': 'vz',
                                 'vr_bc': 'vr', 'mul': 'mu_lcosb', 'mub': 'mu_b',
                                 'pop': 'popid',
                                 'b': 'glat', 'l': 'glon', 'Dist': 'rad',
                                 'log_g': 'grav', 'log_Teff': 'teff', 'Fe/H_initial': 'feh',
                                 'n_companions':'N_companions', 'system_idx':'obj_id', 'system_Mass':'systemMass',
                                 'log_L':'L'},
                         inplace=True)
        companion_df.rename(columns={'iMass':'zams_mass', 'log_Teff': 'teff', 'log_L':'L',
                            'log_g':'logg', 'Mass':'mass', 'Fe/H_initial':'metallicity',
                            'eccentricity':'e'}, inplace=True, errors='ignore')

        system_df['L'] = 10**system_df['L']  # Delogging teff and L
        system_df['teff'] = 10**system_df['teff']
        companion_df['L'] = 10**companion_df['L']
        companion_df['teff'] = 10**companion_df['teff']

        self.logger.info("Basic columns renamed")
        pd.eval('age = log10(system_df.age*1e9)', target=system_df, inplace=True)
        #system_df.loc[:,'age'] = np.log10(system_df['age']*1e9)
        self.logger.info("Age units converted")

        combined_mass = companion_df["zams_mass"] + companion_df["zams_mass_prim"]

        if 'log_a' not in companion_df.columns:
            semimajor_axis = np.cbrt(companion_df["period"]** 2 * combined_mass)
            companion_df["log_a"] = np.log10(semimajor_axis)
        x = np.random.uniform(0, 1, size=companion_df.shape[0])
        y = np.random.uniform(0, 1, size=companion_df.shape[0])
        z = np.random.uniform(0, 1, size=companion_df.shape[0])
        companion_df["i"] = np.arccos(x) * 180 / np.pi
        companion_df["Omega"] = 360 * y
        companion_df["omega"] = 360 * z
        
        #wrap_idx = system_df[system_df['glon'] > 180].index
        #system_df.loc[wrap_idx, 'glon'] -= 360
        pd.eval('glon = system_df.glon - (system_df.glon>180)*360', target=system_df, inplace=True)
        pd.eval("mbol = -2.5 * system_df.L + 4.75", target=system_df, inplace=True)
        # system_df.loc[:, 'mbol2'] = -2.5 * system_df["log_L"].to_numpy() + 4.75
        # self.logger.info("Added mbol2 via loc")
        # system_df.loc[:, 'systemMass2'] = system_df['mass']
        # self.logger.info("Added systemMass2 via loc insertion")
        """pd.eval("systemMass = system_df.mass", target=system_df, inplace=True)
        self.logger.info("Added systemMass via eval")"""

        if system_df['isMultiple'].any():
            map = system_df.set_index('obj_id')['glat'].squeeze()
            companion_df['glat'] = companion_df['system_idx'].map(map)
            map = system_df.set_index('obj_id')['glon'].squeeze()
            companion_df['glon'] = companion_df['system_idx'].map(map)
            
            # system_df.set_index('obj_id')
            # map = companion_df.groupby('system_idx')['mass'].sum()
            # system_df['systemMass'] = system_df['obj_id'].map(map).fillna(0) + system_df['mass']
            # map = companion_df.set_index('system_idx')['mass'].squeeze()
            # system_df['systemMass'] = system_df['obj_id'].map(map)
            # system_df['systemMass'] = system_df['systemMass'] + system_df['mass']
            
        else:
            pd.eval("systemMass = system_df.mass", target=system_df, inplace=True)

        system_df.rename(columns={filter_matching_mist[f]:f for f in self.mag_cols},
                         inplace=True)

        companion_df.rename(columns={filter_matching_mist[f]:f for f in self.mag_cols},
                         inplace=True, errors='ignore')
        
        # system_df.loc[:, 'isMultiple'] = np.zeros(system_df.shape[0], dtype=int)
        # system_df.loc[:, 'N_companions'] = np.zeros(system_df.shape[0], dtype=int)

        phases = np.nan_to_num(system_df['phase'].to_numpy())
        phases[phases == 10] = 101
        phases = phases.astype(int)
        system_df.loc[:, 'rem_id'] = (phases*(phases>100)).astype(int)
        # system_df.loc[:, 'obj_id'] = np.arange(0, len(system_df))

        _, lat_bin_edges, long_bin_edges = _get_bin_edges(latitude, longitude, surveyArea, self.bin_edges_number)

        cols_to_cut = []
        for col in system_df.keys():
            if col not in (popsycle_nonmag_cols + self.mag_cols + list(roman_filters_dict)):
                cols_to_cut.append(col)
        popsycle_df = system_df.drop(columns=cols_to_cut)

        cols_to_cut = []
        for col in companion_df.keys():
            if col not in (popsycle_nonmag_bin_cols + self.mag_cols + list(roman_filters_dict.values())):
                cols_to_cut.append(col)
        popsycle_bin_df = companion_df.drop(columns=cols_to_cut)

        nan_sorted_columns = popsycle_bin_df.isna().sum().sort_values(ascending=False).index  # Sorts column names by amount of NaNs
        popsycle_bin_df = popsycle_bin_df.reindex(columns=nan_sorted_columns)  # Resorts columns to put NaN filled columns at the end
        popsycle_bin_df = popsycle_bin_df.loc[:,~popsycle_bin_df.columns.duplicated()].copy()  # Removes duplicated column names

        # mag cols formatted for PopSyCLE companion dataframe
        filters_ukirt = {"H", "J", "K"}
        filters_ubv = {"U", "B", "V", "R", "I"}
        for filter in filters_ukirt:
            if f"ubv_{filter}" in popsycle_bin_df.columns:
                popsycle_bin_df.rename(columns={f"ubv_{filter}": f"m_ukirt_{filter}"}, inplace=True)
            else:
                popsycle_bin_df[f"m_ukirt_{filter}"] = np.nan

        for filter in filters_ubv:
            if f"ubv_{filter}" in popsycle_bin_df.columns:
                popsycle_bin_df.rename(columns={f"ubv_{filter}": f"m_ubv_{filter}"}, inplace=True)
            else:
                popsycle_bin_df[f"m_ubv_{filter}"] = np.nan

        if self.binning_procedure:
            return popsycle_df, popsycle_bin_df
        else:
            if system_df['isMultiple'].any():
                with h5py.File(f"{self.output_root}_companions.h5", 'w') as h5file:
                    h5file['lat_bin_edges'] = lat_bin_edges
                    h5file['long_bin_edges'] = long_bin_edges
                    
                _bin_lb_hdf5(lat_bin_edges, long_bin_edges, popsycle_bin_df, f"{self.output_root}_companions")
            
            with h5py.File(f"{self.output_root}.h5", 'w') as h5file:
                h5file['lat_bin_edges'] = lat_bin_edges
                h5file['long_bin_edges'] = long_bin_edges
                
            _bin_lb_hdf5(lat_bin_edges, long_bin_edges, popsycle_df, self.output_root)
            self.logger.info(f"PopSyCLE formatted output saved in {self.output_root}.h5")

            if system_df['isMultiple'].any():
                return system_df, companion_df
            else:
                return system_df, None
