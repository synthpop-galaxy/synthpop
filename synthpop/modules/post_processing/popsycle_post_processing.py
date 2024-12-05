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
try:
    from popsycle import ebf
except:
    import ebf

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
        dataframe["Av"] = extinction.extinction_at_lambda(0.544579, ext_in_map)
        R_V = 2.5
        dataframe["exbv"] = R_V * 1/dataframe["Av"]

        ebf_df = pd.DataFrame()
        
        # create log (with same info as galaxia log)
        dtype = [('latitude', 'f8'), ('longitude', 'f8'), ('surveyArea', 'f8')]
        log = np.zeros(1, dtype=dtype)
        log['latitude'] = self.model.l_deg
        log['longitude'] = self.model.b_deg
        log['surveyArea'] = self.model.solid_angle
        
        ebf_df["/popid"] = dataframe["pop"]
        ebf_df["/vx"] = dataframe["U"]
        ebf_df["/vy"] = dataframe["V"]
        ebf_df["/vz"] = dataframe["W"]
        ebf_df["/rad"] = dataframe["Dist"]
        ebf_df["/glat"] = dataframe["b"]
        ebf_df["/glon"] = dataframe["l"]
        ebf_df["/smass"] = dataframe["iMass"]
        ebf_df["/mact"] = dataframe["Mass"]
        ebf_df["/px"] = dataframe["x"]
        ebf_df["/py"] = dataframe["y"]
        ebf_df["/pz"] = dataframe["z"]
        ebf_df["/grav"] = dataframe["logg"]
        ebf_df["/teff"] = dataframe["logTeff"]
        ebf_df["/mbol"] = -2.5 * dataframe["logL"] + 4.75
        ebf_df["/feh"] = dataframe["Fe/H_evolved"]
        ebf_df["/ubv_J"] = dataframe["2MASS_J"]
        ebf_df["/ubv_H"] = dataframe["2MASS_H"]
        ebf_df["/ubv_K"] = dataframe["2MASS_Ks"]
        ebf_df["/ubv_U"] = dataframe["Bessell_U"]
        ebf_df["/ubv_I"] = dataframe["Bessell_I"]
        ebf_df["/ubv_B"] = dataframe["Bessell_B"]
        ebf_df["/ubv_V"] = dataframe["Bessell_V"]
        ebf_df["/ubv_R"] = dataframe["Bessell_R"]
        ebf_df["/exbv_schegel"] = dataframe["E(B-V)"]
        
        # Create a dictionary or structured array with your data
        # Define the output file path
        output_path = self.model.parms.name_for_output + '.ebf'

        # Write data to EBF file
        ebf.write(output_path, "/log", log, 'w')

        for tag in ebf_df:
            try:
                ebf.write(output_path, tag, ebf_df[tag], 'a')
            except:
                print("Error writing", tag)

        return ebf_df
