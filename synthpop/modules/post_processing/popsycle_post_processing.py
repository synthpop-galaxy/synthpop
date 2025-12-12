"""
Post-processing to convert the output into the PopSyCLE input format,
ready to be plugged in at the calc_events stage.
Note: obsmag must be set to FALSE in config file to use this module
"""

__author__ = "M.J. Huston, A. Kim, C. Lee, S. Brooke"

from ._post_processing import PostProcessing
import time
import pandas as pd
import numpy as np
import h5py
import os
from popsycle.synthetic import _get_bin_edges, _bin_lb_hdf5

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

synthpop_nonmag_cols = ['l', 'b', 'Dist',
                        'x', 'y', 'z',
                        'vr_bc', 'mul','mub',
                        'U', 'V', 'W',
                        'iMass','Mass',
                        'log_L', 'log_g', 'log_Teff', '[Fe/H]', 'age',
                        'pop', 'phase', 'n_companions'
                        ]

synthpop_nonmag_bin_cols = ['system_idx', 'period', 'eccentricity', '2MASS_Ks',
                           'star_mass', 'log_R']

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
                        'mbol', 'grav', 'teff', 'feh', 'age',
                        'exbv', 'popid',
                        'isMultiple', 'N_companions', 'rem_id', 'obj_id'] 

class PopsyclePostProcessing(PostProcessing):

    def __init__(self, model, logger, bin_edges_number=None, filter_sets=['ubv'], **kwargs):
        super().__init__(model, logger, **kwargs)
        self.bin_edges_number = bin_edges_number
        #self.filter_sets = filter_sets
        self.mag_cols = []
        self.synthpop_mag_cols = []
        for fset in filter_sets:
            self.mag_cols += [fset+'_'+f for f in filter_set_dict[fset]]
            self.synthpop_mag_cols += [filter_matching_mist[fset+'_'+f] for f in filter_set_dict[fset]]

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Converts DataFrame into format needed for input to PopSyCLE as a replacement
        for Galaxia, saving the file to the set file name + '_psc.h5'.
        """
        self.logger.info(f"Beginning PopSyCLE postprocessing.")
        self.logger.info("DOES NOT ACCOUNT FOR COMPANION STARS CURRENTLY")


        # Drop unused data
        cols_to_cut = []
        for col in system_df.keys():
            if col not in (synthpop_nonmag_cols+self.synthpop_mag_cols+synthpop_nonmag_bin_cols+
                            [self.model.populations[0].extinction.A_or_E_type]):
                cols_to_cut.append(col)
        system_df.drop(columns=cols_to_cut, inplace=True)
        self.logger.info("Unused columns dropped")

        self.output_root = f"{self.model.get_filename(self.model.l_deg, self.model.b_deg)}_psc"

        system_df['isMultiple'] = system_df['n_companions']
        system_df.loc[system_df['isMultiple'] != 0, 'isMultiple'] = 1

        map = system_df.set_index('system_idx')['iMass'].squeeze()
        companion_df['zams_mass_prim'] = companion_df['system_idx'].map(map)

        map = system_df.set_index('system_idx')['b'].squeeze()
        companion_df['glat'] = companion_df['system_idx'].map(map)
        map = system_df.set_index('system_idx')['l'].squeeze()
        companion_df['glon'] = companion_df['system_idx'].map(map)
        
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
            self.logger.info("E(B-V) calculated")
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
        with open(self.output_root + '_synthpop_params.txt', 'w') as params_file:
            params_file.write(f"seed {self.model.parms.random_seed}\n")
        self.logger.info("Parameter file written")

        system_df.rename(columns={'iMass': 'zams_mass','Mass': 'mass',
                                 'x': 'px', 'y': 'py', 'z': 'pz',
                                 'U': 'vx', 'V': 'vy', 'W': 'vz',
                                 'vr_bc': 'vr', 'mul': 'mu_lcosb', 'mub': 'mu_b',
                                 'pop': 'popid',
                                 'b': 'glat', 'l': 'glon', 'Dist': 'rad',
                                 'log_g': 'grav', 'log_Teff': 'teff', '[Fe/H]': 'feh',
                                 'n_companions':'N_companions'}, 
                         inplace=True)
        companion_df.rename(columns={'iMass':'zams_mass', 'log_Teff': 'teff', 'log_L':'L',
                            'log_g':'logg', 'Mass':'mass', '[Fe/H]':'metallicity',
                            'eccentricity':'e'}, inplace=True)

        self.logger.info("Basic columns renamed")
        pd.eval('age = log10(system_df.age*1e9)', target=system_df, inplace=True)
        #system_df.loc[:,'age'] = np.log10(system_df['age']*1e9)
        self.logger.info("Age units converted")

        combined_mass = companion_df["zams_mass"] + companion_df["zams_mass_prim"]
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
        self.logger.info("Galactic longitude wrapped")
        pd.eval("mbol = -2.5 * system_df.log_L + 4.75", target=system_df, inplace=True)
        self.logger.info("Added mbol via eval")
        # system_df.loc[:, 'mbol2'] = -2.5 * system_df["log_L"].to_numpy() + 4.75
        # self.logger.info("Added mbol2 via loc")
        # system_df.loc[:, 'systemMass2'] = system_df['mass']
        # self.logger.info("Added systemMass2 via loc insertion")
        pd.eval("systemMass = system_df.mass", target=system_df, inplace=True)
        self.logger.info("Added systemMass via eval")

        system_df.rename(columns={filter_matching_mist[f]:f for f in self.mag_cols},
                         inplace=True)

        companion_df.rename(columns={filter_matching_mist[f]:f for f in self.mag_cols},
                         inplace=True)
        
        self.logger.info("Renamed mag cols")
        
        # system_df.loc[:, 'isMultiple'] = np.zeros(system_df.shape[0], dtype=int)
        # system_df.loc[:, 'N_companions'] = np.zeros(system_df.shape[0], dtype=int)

        phases = np.nan_to_num(system_df['phase'].to_numpy())
        system_df.loc[:, 'rem_id'] = (phases*(phases>100)).astype(int)
        system_df.loc[:, 'obj_id'] = np.arange(0, len(system_df))
        self.logger.info("Multiplicity and remnant columns added")

        self.logger.info("PopSyCLE system_df modifications complete - binning and saving file")
        _, lat_bin_edges, long_bin_edges = _get_bin_edges(latitude, longitude, surveyArea, self.bin_edges_number)

        cols_to_cut = []
        for col in system_df.keys():
            if col not in (popsycle_nonmag_cols + self.mag_cols):
                cols_to_cut.append(col)
        popsycle_df = system_df.drop(columns=cols_to_cut)

        cols_to_cut = []
        for col in companion_df.keys():
            if col not in (popsycle_nonmag_bin_cols + self.mag_cols):
                cols_to_cut.append(col)
        popsycle_bin_df = companion_df.drop(columns=cols_to_cut)

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

        with h5py.File(f"{self.output_root}_companions.h5", 'w') as h5file:
            h5file['lat_bin_edges'] = lat_bin_edges
            h5file['long_bin_edges'] = long_bin_edges

        _bin_lb_hdf5(lat_bin_edges, long_bin_edges, popsycle_bin_df, f"{self.output_root}_companions")
        
        with h5py.File(f"{self.output_root}.h5", 'w') as h5file:
            h5file['lat_bin_edges'] = lat_bin_edges
            h5file['long_bin_edges'] = long_bin_edges

        _bin_lb_hdf5(lat_bin_edges, long_bin_edges, popsycle_df, self.output_root)
        self.logger.info(f"PopSyCLE formatted output saved in {self.output_root}.h5")
    
        return system_df, companion_df
