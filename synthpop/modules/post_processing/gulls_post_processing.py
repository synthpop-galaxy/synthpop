"""
Post-processing to put catalog in the Gulls input format.
"""

__all__ = ["GullsPostProcessing", ]
__author__ = "Farzaneh Zohrabi, Ali Crisp"
__date__ = "2025-05-15"

from ._post_processing import PostProcessing
import time
import pandas as pd
import numpy as np
from synthpop.synthpop_utils.coordinates_transformation import lb_to_ad

class GullsPostProcessing(PostProcessing):
    def __init__(self, model, cat_type, **kwargs):
        """
        Parameters:
            model: SynthPop
                to access properties from the model.
                ege model.const, model.params, model.populations[0]
            cat_type:
                settings for the catalog type
                must be one of : lens, source, extra-bright, bright, mid1, mid2, faint

            kwargs: dict
                keyword arguments from the config file
        """
        self.cat_type = cat_type
        self.model = model

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):

        if companion_df is not None:
            raise ValueError("Must run combine_tables postproc before"+ \
                                "gulls_post_processing")

        # convert l, b to ra, dec
        system_df["l(deg)"] = system_df.iloc[:, 7]  # l
        system_df["b(deg)"] = system_df.iloc[:, 8]  # b
        ra, dec = lb_to_ad(system_df["l(deg)"], system_df["b(deg)"])
        system_df["RA2000.0"] = ra
        system_df["DEC2000.0"] = dec

        # convert columns
        system_df["logg"] = system_df["log_g"]
        system_df["Mbol"] = -2.5 * system_df["log_L"] + 4.75
        system_df["Teff"] = 10 ** system_df["log_Teff"]
        system_df["[alpha/Fe]"] = 0
        system_df["Radius"] = 10 ** system_df["log_R"]
        system_df["CL"] = system_df["phase"]
        system_df["Vr"] = system_df['vr_bc']
        
        # compute an approximate magnitude for Roman F213 from 2MASS Ks
        k213 = system_df["2MASS_Ks"] + 1.834505
       
        # reduce data frame to the needed columns
        filtlist = self.model.parms.chosen_bands
        cols = filtlist + ["mul", "mub", "Vr", "U", "V", "W", "iMass", 
                "CL", "age", "Teff", "logg", "pop", "Mass", "Mbol", "Radius", 
                "[Fe/H]","l", "b", "RA2000.0", "DEC2000.0", 
                "Dist", "x", "y", "z", "A_Ks", "[alpha/Fe]"]
        system_df = system_df[cols]
        
        # Get index of F184 band and insert K213 after that
        idx = system_df.columns.get_loc("F184")+1
        system_df.insert(idx, "K213", k213)
        
        
        # Replace NaNs with 99 for magnitude columns, and a non-physical
        # very small value for other parameters.
        for filt in filtlist:
            system_df.loc[:, filt] = system_df.loc[:,filt].fillna(99)
        system_df.loc[:, "K213"] = system_df.loc[:,"K213"].fillna(99)
        system_df.loc[:, "Mbol"] = system_df.loc[:,"Mbol"].fillna(99)
        system_df = system_df.fillna(value=2e-50)


        # Impose magnitude lower limits based on Roman expectations
        # Lens has no limits, source and extra bright only have faint limits
        # All faint limits should be handled by the configuration file, so are
        # not set here.
        if self.cat_type == "lens":
            pass

        elif self.cat_type == "source":
            pass
            
        elif self.cat_type == "extra-bright":
            pass
            
        elif self.cat_type == "bright":
            system_df = system_df.loc[system_df["W146"] > 12]

        elif self.cat_type == "mid1":
            system_df = system_df.loc[system_df["W146"] > 17]

        elif self.cat_type == "mid2":
            system_df = system_df.loc[system_df["W146"] > 20]

        elif self.cat_type == "faint":
            system_df = system_df.loc[system_df["W146"] > 24]

        return system_df, companion_df
