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

filter_matching_mist = {"2MASS_J": 'ubv_J',
                   "2MASS_H": 'ubv_H',
                   "2MASS_Ks": 'ubv_K',
                   "Bessell_U": 'ubv_U',
                   "Bessell_I": 'ubv_I',
                   "Bessell_B": 'ubv_B',
                   "Bessell_V": 'ubv_V',
                   "Bessell_R": 'ubv_R'
                   }

class PopsyclePostProcessing(PostProcessing):

    def __init__(self, model, logger, bin_edges_number=None, **kwargs):
        super().__init__(model, logger, **kwargs)
        self.bin_edges_number = bin_edges_number

    def do_post_processing(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Converts DataFrame into format needed for input to PopSyCLE as a replacement
        for Galaxia, saving the file to the set file name + '.ebf' and returning
        the DataFrame with modified columns.
        """
        self.output_root = f"{self.model.get_filename(self.model.l_deg, self.model.b_deg, self.model.solid_angle)}_psc"
        
        # Translate extinction to Ebv extinction
        extinction = self.model.populations[0].extinction
        ext_in_map = dataframe.iloc[:, 19]
        Av = extinction.extinction_at_lambda(0.544579, ext_in_map)
        Ab = extinction.extinction_at_lambda(0.438074, ext_in_map)
        dataframe["E(B-V)"] = Ab - Av

        star_dict = pd.DataFrame()
        
        # create log (with same info as galaxia log)
        #dtype = [('latitude', 'f8'), ('longitude', 'f8'), ('surveyArea', 'f8')]
        #log = np.zeros(1, dtype=dtype)
        latitude = self.model.l_deg
        longitude = self.model.b_deg
        surveyArea = self.model.solid_angle
        if self.model.solid_angle_unit=='sr':
            surveyArea *= (180/np.pi)**2
            
        with open(self.output_root + '_synthpop_params.txt', 'w') as params_file:
            params_file.write(f"seed {self.model.parms.random_seed}\n")

        star_dict['zams_mass'] = dataframe["iMass"]
        star_dict['mass'] = dataframe["Mass"]
        star_dict['systemMass'] = star_dict['mass']
        star_dict['px'] = dataframe["x"]
        star_dict['py'] = dataframe["y"]
        star_dict['pz'] = dataframe["z"]
        star_dict['vx'] = dataframe["U"]
        star_dict['vy'] = dataframe["V"]
        star_dict['vz'] = dataframe["W"]
        star_dict['vr'] = dataframe['vr_bc']
        star_dict['mu_lcosb'] = dataframe['mul']
        star_dict['mu_b'] = dataframe['mub']
        star_dict['exbv'] = dataframe["E(B-V)"]
        star_dict['glat'] = dataframe["b"]
        star_dict['glon'] = dataframe["l"]
        wrap_idx = np.where(star_dict['glon'] > 180)[0]
        star_dict['glon'][wrap_idx] -= 360
        
        star_dict['mbol'] = -2.5 * dataframe["logL"] + 4.75
        star_dict['grav'] = dataframe["logg"]
        star_dict['teff'] = dataframe["logTeff"]
        star_dict['feh'] = dataframe["Fe/H_evolved"]
        star_dict['rad'] = dataframe["Dist"]
        for filter in self.model.parms.chosen_bands:
            if filter in filter_matching_mist:
                star_dict[filter_matching_mist[filter]] = dataframe[filter]
        
        star_dict['isMultiple'] = np.zeros(dataframe.shape[0], dtype=int)
        star_dict['N_companions'] = np.zeros(dataframe.shape[0], dtype=int)
        phases = np.nan_to_num(dataframe['phase'].to_numpy())
        star_dict['rem_id'] = (phases*(phases>100)).astype(int)
        star_dict['obj_id'] = np.arange(0, len(star_dict))

        _, lat_bin_edges, long_bin_edges = _get_bin_edges(latitude, longitude, surveyArea, self.bin_edges_number)
        
        h5file = h5py.File(f"{self.output_root}.h5", 'w')
        h5file['lat_bin_edges'] = lat_bin_edges
        h5file['long_bin_edges'] = long_bin_edges
        
        _bin_lb_hdf5(lat_bin_edges, long_bin_edges, star_dict, self.output_root)
        print(f"PopSyCLE formatted output saved in {self.output_root}.h5")
    
        return dataframe
