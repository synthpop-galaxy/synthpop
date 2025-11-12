"""
Post-processing to convert the output into the PopSyCLE input format,
ready to be plugged in at the calc_events stage.
Note: obsmag must be set to FALSE in config file to use this module
"""

__author__ = "M.J. Huston, A. Kim, C. Lee"

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
                        'pop', 'phase'
                        ]
synthpop_nonmag_bin_cols = []

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

    def do_post_processing(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Converts DataFrame into format needed for input to PopSyCLE as a replacement
        for Galaxia, saving the file to the set file name + '_psc.h5'.
        """
        self.logger.info(f"Beginning PopSyCLE postprocessing.")

        # Drop unused data
        cols_to_cut = []
        for col in dataframe.keys():
            if col not in (synthpop_nonmag_cols+self.synthpop_mag_cols+synthpop_nonmag_bin_cols+
                            [self.model.populations[0].extinction.A_or_E_type]):
                cols_to_cut.append(col)
        dataframe.drop(columns=cols_to_cut, inplace=True)
        self.logger.info("Unused columns dropped")

        self.output_root = f"{self.model.get_filename(self.model.l_deg, self.model.b_deg, self.model.solid_angle)}_psc"
        
        # Translate extinction to Ebv extinction
        if not self.model.populations[0].extinction.A_or_E_type=="E(B-V)":
            extinction = self.model.populations[0].extinction
            extinction_type = self.model.populations[0].extinction.A_or_E_type
            dataframe.rename(columns={extinction_type:'ext_orig'}, inplace=True)
            #ext_in_map = dataframe[extinction_type].to_numpy()
            Av_Aref = extinction.extinction_at_lambda(0.544579, 1.0)
            Ab_Aref = extinction.extinction_at_lambda(0.438074, 1.0)
            pd.eval("exbv = (Ab_Aref - Av_Aref)*dataframe.ext_orig", target=dataframe, inplace=True)
            dataframe.drop(columns=['ext_orig'], inplace=True)
            self.logger.info("E(B-V) calculated")
        else:
            dataframe.rename(columns={"E(B-V)":'exbv'}, inplace=True)
        
        # create log (with same info as galaxia log)
        #dtype = [('latitude', 'f8'), ('longitude', 'f8'), ('surveyArea', 'f8')]
        #log = np.zeros(1, dtype=dtype)
        latitude = self.model.l_deg
        longitude = self.model.b_deg
        surveyArea = self.model.solid_angle
        if self.model.solid_angle_unit=='sr':
            surveyArea *= (180/np.pi)**2
            
        os.makedirs('/'.join(self.output_root.split('/')[:-1]), exist_ok=True)
        with open(self.output_root + '_synthpop_params.txt', 'w') as params_file:
            params_file.write(f"seed {self.model.parms.random_seed}\n")
        self.logger.info("Parameter file written")

        dataframe.rename(columns={'iMass': 'zams_mass','Mass': 'mass', 
                                 'x': 'px', 'y': 'py', 'z': 'pz',
                                 'U': 'vx', 'V': 'vy', 'W': 'vz',
                                 'vr_bc': 'vr', 'mul': 'mu_lcosb', 'mub': 'mu_b',
                                 'pop': 'popid',
                                 'b': 'glat', 'l': 'glon', 'Dist': 'rad',
                                 'log_g': 'grav', 'log_Teff': 'teff', '[Fe/H]': 'feh'}, 
                         inplace=True)
        self.logger.info("Basic columns renamed")
        pd.eval('age = log10(dataframe.age*1e9)', target=dataframe, inplace=True)
        #dataframe.loc[:,'age'] = np.log10(dataframe['age']*1e9)
        self.logger.info("Age units converted")

        #wrap_idx = dataframe[dataframe['glon'] > 180].index
        #dataframe.loc[wrap_idx, 'glon'] -= 360
        pd.eval('glon = dataframe.glon - (dataframe.glon>180)*360', target=dataframe, inplace=True)
        self.logger.info("Galactic longitude wrapped")
        pd.eval("mbol = -2.5 * dataframe.log_L + 4.75", target=dataframe, inplace=True)
        self.logger.info("Added mbol via eval")
        # dataframe.loc[:, 'mbol2'] = -2.5 * dataframe["log_L"].to_numpy() + 4.75
        # self.logger.info("Added mbol2 via loc")
        # dataframe.loc[:, 'systemMass2'] = dataframe['mass']
        # self.logger.info("Added systemMass2 via loc insertion")
        pd.eval("systemMass = dataframe.mass", target=dataframe, inplace=True)
        self.logger.info("Added systemMass via eval")


        dataframe.rename(columns={filter_matching_mist[f]:f for f in self.mag_cols},
                         inplace=True)
        self.logger.info("Renamed mag cols")
        
        dataframe.loc[:, 'isMultiple'] = np.zeros(dataframe.shape[0], dtype=int)
        dataframe.loc[:, 'N_companions'] = np.zeros(dataframe.shape[0], dtype=int)
        phases = np.nan_to_num(dataframe['phase'].to_numpy())
        dataframe.loc[:, 'rem_id'] = (phases*(phases>100)).astype(int)
        dataframe.loc[:, 'obj_id'] = np.arange(0, len(dataframe))
        self.logger.info("Multiplicity and remnant columns added")

        self.logger.info("PopSyCLE dataframe modifications complete - binning and saving file")
        _, lat_bin_edges, long_bin_edges = _get_bin_edges(latitude, longitude, surveyArea, self.bin_edges_number)
        
        with h5py.File(f"{self.output_root}.h5", 'w') as h5file:
            h5file['lat_bin_edges'] = lat_bin_edges
            h5file['long_bin_edges'] = long_bin_edges
        
        _bin_lb_hdf5(lat_bin_edges, long_bin_edges, dataframe[popsycle_nonmag_cols+self.mag_cols], self.output_root)
        self.logger.info(f"PopSyCLE formatted output saved in {self.output_root}.h5")
    
        return dataframe