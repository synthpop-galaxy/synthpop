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
import os

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
        star_dict = pd.DataFrame()
        
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
        star_dict['ubv_J'] = dataframe["VISTA_J"]
        star_dict['ubv_H'] = dataframe["VISTA_H"]
        star_dict['ubv_K'] = dataframe["VISTA_Ks"]
        star_dict['ubv_U'] = dataframe["Bessell_U"]
        star_dict['ubv_I'] = dataframe["Bessell_I"]
        star_dict['ubv_B'] = dataframe["Bessell_B"]
        star_dict['ubv_V'] = dataframe["Bessell_V"]
        star_dict['ubv_R'] = dataframe["Bessell_R"]
        
        star_dict['isMultiple'] = np.zeros(dataframe.shape[0], dtype=int)
        star_dict['N_companions'] = np.zeros(dataframe.shape[0], dtype=int)
        star_dict['rem_id'] = dataframe["Dim_Compact_Object_Flag"].map({0.0: 0, 1.0: 101, 2.0: 102, 3.0: 103})
        
#         print("px range:", star_dict["px"].min(), star_dict["px"].max())
#         print("vx range:", star_dict["vx"].min(), star_dict["vx"].max())
#         print("mass range:", star_dict["mass"].min(), star_dict["mass"].max())

        
        # Using only 1 bin
        l_min, l_max = dataframe["l"].min(), dataframe["l"].max()
        b_min, b_max = dataframe["b"].min(), dataframe["b"].max()
        
        delta_l = 0.1 if l_min == l_max else 0
        delta_b = 0.1 if b_min == b_max else 0

        long_bin_edges = np.array([l_min, l_max + delta_l])
        lat_bin_edges = np.array([b_min, b_max + delta_b])

#         long_bin_edges = np.array([l_min, l_max])
#         lat_bin_edges = np.array([b_min, b_max])
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
    
        # Using h5 file creation code from PopSyCLE
        def lb_binning_hdf5(obj_arr: pd.DataFrame, output_root: str, n_lat_bins=2, n_lon_bins=2, companion_obj_arr=None):

            file_path = output_root + '.h5'
        
            # Determine lat-lon ranges and edges
            l_min, l_max = obj_arr['glon'].min(), obj_arr['glon'].max()
            b_min, b_max = obj_arr['glat'].min(), obj_arr['glat'].max()
        
            long_bin_edges = np.linspace(l_min, l_max, n_lon_bins + 1)
            lat_bin_edges = np.linspace(b_min, b_max, n_lat_bins + 1)
        
            # Convert to structured array and dtype
            compound_dtype = np.dtype([(col, obj_arr[col].dtype) for col in obj_arr.columns])
            obj_arr = obj_arr.to_records(index=False)
        
            if companion_obj_arr is not None:
                compound_dtype = np.dtype(companion_obj_arr)
        
            # Create file and store bin edges
            if not os.path.exists(file_path):
                with h5py.File(file_path, 'w') as hf:
                    hf.create_dataset("lat_bin_edges", data=lat_bin_edges)
                    hf.create_dataset("long_bin_edges", data=long_bin_edges)
        
            # Bin loop
            for ll in range(n_lon_bins):
                for bb in range(n_lat_bins):
                    with h5py.File(file_path, 'r+') as hf:
                        dset_name = f"l{ll}b{bb}"
        
                        # Filter for lat-long bin
                        if companion_obj_arr is None:
                            in_bin = obj_arr[
                                (obj_arr['glat'] >= lat_bin_edges[bb]) & (obj_arr['glat'] < lat_bin_edges[bb + 1]) &
                                (obj_arr['glon'] >= long_bin_edges[ll]) & (obj_arr['glon'] < long_bin_edges[ll + 1])
                            ]
        
                            if len(in_bin) == 0:
                                continue
        
                            save_data = np.array(in_bin, dtype=compound_dtype)
        
                        else:
                            id_lb = np.where(
                                (obj_arr['glat'] >= lat_bin_edges[bb]) & (obj_arr['glat'] < lat_bin_edges[bb + 1]) &
                                (obj_arr['glon'] >= long_bin_edges[ll]) & (obj_arr['glon'] < long_bin_edges[ll + 1])
                            )[0]
        
                            if len(id_lb) == 0:
                                continue
        
                            primary_ids = obj_arr['obj_id'][id_lb]
                            companion_id_lb = [np.where(companion_obj_arr['system_idx'] == pid)[0] for pid in primary_ids]
                            companion_id_lb = np.concatenate(companion_id_lb)
        
                            if len(companion_id_lb) == 0:
                                continue
        
                            save_data = np.array(companion_obj_arr[companion_id_lb], dtype=compound_dtype)
        
                        # Create or append to dataset
                        if dset_name not in hf:
                            dset = hf.create_dataset(
                                dset_name,
                                data=save_data,
                                maxshape=(None,),
                                chunks=(min(10000, len(save_data)),),
                                dtype=compound_dtype,
                                compression="gzip"
                            )
                        else:
                            dset = hf[dset_name]
                            old_size = dset.shape[0]
                            new_size = old_size + len(save_data)
                            dset.resize((new_size,))
                            dset[old_size:new_size] = save_data
        
            return file_path

        
        lb_binning_hdf5(lat_bin_edges, long_bin_edges, star_dict, "synthpop_test_h5_file")
    
        return dataframe
