"""
NOTE: IN PROGRESS MODULE - NOT READY FOR USE
Post-processing to convert the output into the PopSyCLE input format,
as a replacement for a Galaxia catalog.
Note: obsmag must be set to FALSE in config file to use this module
"""

__author__ = "A. Kim, M.J. Huston"

from ._post_processing import PostProcessing
import time
import pandas as pd
import numpy as np
import h5py

class PopsyclePostProcessing(PostProcessing):

    def __init__(self, model, logger, **kwargs):
        super().__init__(model, logger, **kwargs)

    def do_post_processing(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Converts DataFrame into format needed for input to PopSyCLE as a replacement
        for Galaxia, saving the file to the set file name + '.ebf' and returning
        the DataFrame with modified columns.
        """
        # Translate extinction to Ebv extinction
        extinction = self.model.populations[0].extinction
        ext_in_map = dataframe.iloc[:, 19]
        Av = extinction.extinction_at_lambda(0.544579, ext_in_map)
        Ab = extinction.extinction_at_lambda(0.438074, ext_in_map)
        dataframe["E(B-V)"] = Ab - Av

#        ebf_df = pd.DataFrame()
        
        # create log (with same info as galaxia log)
        dtype = [('latitude', 'f8'), ('longitude', 'f8'), ('surveyArea', 'f8')]
        log = np.zeros(1, dtype=dtype)
        log['latitude'] = self.model.l_deg
        log['longitude'] = self.model.b_deg
        log['surveyArea'] = self.model.solid_angle
        
#        ebf_df["/popid"] = dataframe["pop"]
#        ebf_df["/vx"] = dataframe["U"]
#        ebf_df["/vy"] = dataframe["V"]
#        ebf_df["/vz"] = dataframe["W"]
#        ebf_df["/rad"] = dataframe["Dist"]
#        ebf_df["/glat"] = dataframe["b"]
#        ebf_df["/glon"] = dataframe["l"]
#        ebf_df["/smass"] = dataframe["iMass"]
#        ebf_df["/mact"] = dataframe["Mass"]
#        ebf_df["/px"] = dataframe["x"]
#        ebf_df["/py"] = dataframe["y"]
#        ebf_df["/pz"] = dataframe["z"]
#        ebf_df["/grav"] = dataframe["logg"]
#        ebf_df["/teff"] = dataframe["logTeff"]
#        ebf_df["/mbol"] = -2.5 * dataframe["logL"] + 4.75
#        ebf_df["/feh"] = dataframe["Fe/H_evolved"]
#        ebf_df["/ubv_J"] = dataframe["2MASS_J"]
#        ebf_df["/ubv_H"] = dataframe["2MASS_H"]
#        ebf_df["/ubv_K"] = dataframe["2MASS_Ks"]
#        ebf_df["/ubv_U"] = dataframe["Bessell_U"]
#        ebf_df["/ubv_I"] = dataframe["Bessell_I"]
#        ebf_df["/ubv_B"] = dataframe["Bessell_B"]
#        ebf_df["/ubv_V"] = dataframe["Bessell_V"]
#        ebf_df["/ubv_R"] = dataframe["Bessell_R"]
#        ebf_df["/exbv_schegel"] = dataframe["E(B-V)"]

        star_dict['zams_mass'] = dataframe["iMass"]
        star_dict['mass'] = dataframe["Mass"]
        star_dict['systemMass'] = star_dict['mass']
        star_dict['px'] = dataframe["x"]
        star_dict['py'] = dataframe["y"]
        star_dict['pz'] = dataframe["z"]
        star_dict['vx'] = dataframe["U"]
        star_dict['vy'] = dataframe["V"]
        star_dict['vz'] = dataframe["W"]
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
        star_dict['ubv_J'] = dataframe["2MASS_J"]
        star_dict['ubv_H'] = dataframe["2MASS_H"]
        star_dict['ubv_K'] = dataframe["2MASS_Ks"]
        star_dict['ubv_U'] = dataframe["Bessell_U"]
        star_dict['ubv_I'] = dataframe["Bessell_I"]
        star_dict['ubv_B'] = dataframe["Bessell_B"]
        star_dict['ubv_V'] = dataframe["Bessell_V"]
        star_dict['ubv_R'] = dataframe["Bessell_R"]
        
        star_dict['isMultiple'] = np.zeros(len(dataframe.shape[0]), dtype=int)
        star_dict['N_companions'] = np.zeros(len(dataframe.shape[0]), dtype=int)
        star_dict['rem_id'] = dataframe["Dim_Compact_Object_Flag"].map({0.0: 0, 1.0: 101, 2.0: 102, 3.0: 103})
        
        
        # Using only 1 bin
        l_min, l_max = dataframe["l"].min(), dataframe["l"].max()
        b_min, b_max = dataframe["b"].min(), dataframe["b"].max()
        long_bin_edges = np.array([l_min, l_max])
        lat_bin_edges = np.array([b_min, b_max])
        bin_edges_number = 2
        
        
        #         # Create a dictionary or structured array with your data
#         # Define the output file path
#         output_path = self.model.parms.name_for_output + '.ebf'

#         # Write data to EBF file
#         ebf.write(output_path, "/log", log, 'w')

#         for tag in ebf_df:
#             try:
#                 ebf.write(output_path, tag, ebf_df[tag], 'a')
#             except:
#                 print("Error writing", tag)

#         return ebf_df

        lb_binning_hdf5(lat_bin_edges, long_bin_edges, star_dict, "synthpop_test_h5_file")
    
        # Using h5 file creation code from PopSyCLE
        def lb_binning_hdf5(lat_bin_edges, long_bin_edges, obj_arr, output_root, companion_obj_arr = None):

    
            if companion_obj_arr is None:
                # Create compound datatype from obj_arr
                compound_dtype = _generate_compound_dtype(obj_arr)
            else:
                compound_dtype = np.dtype(companion_obj_arr)

            ##########
            # Loop through the latitude and longitude bins.
            ##########
            
            for ll in range(len(long_bin_edges) - 1):
                for bb in range(len(lat_bin_edges) - 1):
                    # Open our HDF5 file for reading and appending.
                    # Create as necessary.
                    hf = h5py.File(output_root + '.h5', 'r+')

                    # HDF5 dataset name
                    dset_name = 'l' + str(ll) + 'b' + str(bb)

                    # Create data set if needed. Start with 0 stars in the dataset.
                    if dset_name not in hf:
                        dataset = hf.create_dataset(dset_name, shape=(0,),
                                                    chunks=(1e4,),
                                                    maxshape=(None,),
                                                    dtype=compound_dtype)
                    else:
                        dataset = hf[dset_name]

                    ##########
                    # Binning the stars and/or compact objects or companions
                    ##########
                    if obj_arr is not None:
                        id_lb = np.where((obj_arr['glat'] >= lat_bin_edges[bb]) &
                                         (obj_arr['glat'] < lat_bin_edges[bb + 1]) &
                                         (obj_arr['glon'] >= long_bin_edges[ll]) &
                                         (obj_arr['glon'] < long_bin_edges[ll + 1]))[0]

                        if len(id_lb) == 0:
                            continue

                        # Loop over the obj_arr and add all columns
                        # (matching id_lb) into save_data
                        save_data = np.empty(len(id_lb), dtype=compound_dtype)

                        if companion_obj_arr is None:
                            for colname in obj_arr:
                                save_data[colname] = obj_arr[colname][id_lb]
                        # If making a companion hd5f file, finds corresponding companions and save them
                        else:
                            companion_id_lb = [np.where(companion_obj_arr['system_idx'] == ii)[0] for ii in obj_arr['obj_id'][id_lb]]
                            companion_id_lb = list(np.concatenate(companion_id_lb).ravel()) # Simplifies datastructure
                            if len(companion_id_lb) == 0:
                                continue
                            save_data = np.array(companion_obj_arr[companion_id_lb])

                        # Resize the dataset and add data.
                        old_size = dataset.shape[0]
                        if companion_obj_arr is None:
                            new_size = old_size + len(id_lb)                
                        else:
                            new_size = old_size + len(companion_id_lb)                     
                        dataset.resize((new_size, ))
                        dataset[old_size:new_size] = save_data

                    hf.close()
                    
            return output_root.h5
