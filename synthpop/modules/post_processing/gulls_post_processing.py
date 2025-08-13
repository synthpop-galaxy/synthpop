"""
Post-processing to put catalog in the Gulls input format.
"""

__all__ = ["GullsPostProcessing", ]
__author__ = "Farzaneh Zohrabi, Ali Crisp"
__date__ = "2025-05-15"

from ._post_processing import PostProcessing
import time
import pandas
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

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """

        Parameters
        ----------
        dataframe : dataframe
            original SynthPop output as pandas data frame

        Returns
        -------
        dataframe : dataframe
            modified pandas data frame
            
        """

        # convert l, b to ra, dec
        dataframe["l(deg)"] = dataframe.iloc[:, 7]  # l
        dataframe["b(deg)"] = dataframe.iloc[:, 8]  # b
        ra, dec = lb_to_ad(dataframe["l(deg)"], dataframe["b(deg)"])
        dataframe["RA2000.0"] = ra
        dataframe["DEC2000.0"] = dec

        # convert columns
        dataframe["logg"] = dataframe["log_g"]
        dataframe["Mbol"] = -2.5 * dataframe["log_L"] + 4.75
        dataframe["Teff"] = 10 ** dataframe["log_Teff"]
        dataframe["[alpha/Fe]"] = 0
        dataframe["Radius"] = 10 ** dataframe["log_R"]
        dataframe["CL"] = dataframe["phase"]
        dataframe["Vr"] = dataframe['vr_bc']
        
        # compute an approximate magnitude for Roman F213 from 2MASS Ks
        k213 = dataframe["2MASS_Ks"] + 1.834505
       
        # reduce data frame to the needed columns
        filtlist = self.model.parms.chosen_bands
        cols = filtlist + ["mul", "mub", "Vr", "U", "V", "W", "iMass", 
                "CL", "age", "Teff", "logg", "pop", "Mass", "Mbol", "Radius", 
                "[Fe/H]","l", "b", "RA2000.0", "DEC2000.0", 
                "Dist", "x", "y", "z", "A_Ks", "[alpha/Fe]"]
        dataframe = dataframe[cols]
        
        # Get index of F184 band and insert K213 after that
        idx = dataframe.columns.get_loc("F184")+1
        dataframe.insert(idx, "K213", k213)
        
        
        # Replace NaNs with 99 for magnitude columns, and a non-physical
        # very small value for other parameters.
        for filt in filtlist:
            dataframe.loc[:, filt] = dataframe.loc[:,filt].fillna(99)
        dataframe.loc[:, "K213"] = dataframe.loc[:,"K213"].fillna(99)
        dataframe.loc[:, "Mbol"] = dataframe.loc[:,"Mbol"].fillna(99)
        dataframe = dataframe.fillna(value=2e-50)


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
            dataframe = dataframe.loc[dataframe["W146"] > 12]

        elif self.cat_type == "mid1":
            dataframe = dataframe.loc[dataframe["W146"] > 17]

        elif self.cat_type == "mid2":
            dataframe = dataframe.loc[dataframe["W146"] > 20]

        elif self.cat_type == "faint":
            dataframe = dataframe.loc[dataframe["W146"] > 24]

        return dataframe
