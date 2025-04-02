"""
Post-processing to put catalog in the Gulls input format.
"""

__all__ = ["GullsPostProcessing", ]
__author__ = "Farzaneh Zohrabi"

from ._post_processing import PostProcessing
import time
import pandas
import numpy as np
from synthpop.synthpop_utils.coordinates_transformation import lb_to_ad

class GullsPostProcessing(PostProcessing):
    @staticmethod
    def Alambda_AV(eff_wavelength: float, R_V: float = 2.5) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the relative extinction A_lambda/A_V
    
        Parameters
        ----------
        eff_wavelength : float
            Effective Wavelength of the filter for which the extinction should be determined.
            in micrometer
        R_V : float
            interstellar reddening parameter
        """
    
        x = 1 / eff_wavelength
        # deep red?
        # if x <0.3:
        #    return (1/x)**(-1.75)
    
        # infrared
        # if x>= 0.3 and x<=1.1:
        if x <= 1.1:
            a = 0.574 * x ** 1.61
            b = -0.527 * x ** 1.61
    
        # optical
        elif 1.1 < x <= 3.3:
            y = x - 1.82
            a = (1 + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4
                 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7)
            b = (1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4
                 - 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7)
    
        # UV and far UV
        elif 3.3 < x <= 8:
            if x < 5.9:
                F_a = 0
                F_b = 0
            else:
                F_a = -0.04473 * (x - 5.9) ** 2 - 0.009779 * (x - 5.9) ** 3.
                F_b = 0.2130 * (x - 5.9) ** 2 + 0.1207 * (x - 5.9) ** 3.
            a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) + F_a
            b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) + F_b
        else:
            a, b = None, None
    
        return a + b / R_V
    @staticmethod
    def Alambda_Afilt(eff_wavelength: float, R_V: float = 2.5):
        """
        Calculate the extinction relative to the specified filter
        for a given effective wavelength.
                 
        Arguments
        ---------
        eff_wavelength : float
            Effective Wavelength of the filter for which the extinction should be determined.
            in micrometer
        R_V : float
            interstellar reddening parameter
        """
        ref_wavelength = 2.152152
        Afilt_AV = GullsPostProcessing.Alambda_AV(ref_wavelength, R_V)
        AlambdaAV = GullsPostProcessing.Alambda_AV(eff_wavelength, R_V)
        # Return  A_lambda   A_lambda    A_V
        #        -------- = -------- * ------
        #         A_filt      A_V      A_filt
        return AlambdaAV / Afilt_AV
    
    def __init__(self, model, cat_type, **kwargs):
        """
        Parameters:
            model: SynthPop
                to access properties from the model.
                ege model.const, model.params, model.populations[0]
            cat_type:
                settings for the catalog type
                must be on of :

            kwargs: dict
                keyword arguments from the config file
        """
        self.cat_type = cat_type
        self.model = model

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        This is a placeholder for the postprocessing
        replace it whit what ever you think is useful
        It must accept the pandas data frame and must return a pandas data frame

        Parameters
        ----------
        dataframe : dataframe
            original SynthPop output as pandas data frame

        Returns
        -------
        dataframe : dataframe
            modified pandas data frame
            
        """
        # Translate extinction to V band extinction
        extinction = self.model.populations[0].extinction

        ext_in_map = dataframe.iloc[:, 19]
        Alambda_A_map = self.Alambda_Afilt(0.544579,2.5)
        dataframe["Av"]= dataframe["A_Ks"]* Alambda_A_map
        # convert l, b to ra, dec
        dataframe["l(deg)"] = dataframe.iloc[:, 7]  # l
        dataframe["b(deg)"] = dataframe.iloc[:, 8]  # b
        ra, dec = lb_to_ad(dataframe["l(deg)"], dataframe["b(deg)"])
        dataframe["RA2000.0"] = ra
        dataframe["DEC2000.0"] = dec

        # convert columns
        dataframe["[Fe/H]"] = dataframe["Fe/H_evolved"]
        dataframe["Mbol"] = -2.5 * dataframe["logL"] + 4.75
        dataframe["Teff"] = 10 ** dataframe["logTeff"]
        dataframe["[alpha/Fe]"] = 0
        dataframe["Radius"] = 10 ** dataframe["log_radius"]
        #dataframe["Cl"] = dataframe["phase"]
        dataframe["Vr"]=dataframe['vr_bc']
        # converting Vega mag to AB mag
        # Vega_AB_Correction = {
        #     "R062": 0.137095,
        #     "Z087": 0.487379,
        #     "Y106": 0.653780,
        #     "J129": 0.958363,
        #     "W146": 1.024467,
        #     "H158": 1.287404,
        #     "F184": 1.551332,
        #     "2MASS_Ks": 1.834505,
        #     "2MASS_J": 0.889176,
        #     "2MASS_H": 1.364157,
        #     "Gaia_G_EDR3": 0.128961,
        #     "Gaia_BP_EDR3": 0.030233,
        #     "Gaia_RP_EDR3": 0.373349,
        #     "Bessell_U": 0.800527,
        #     "Bessell_B": -0.107512,
        #     "Bessell_V": 0.006521,
        #     "Bessell_I": 0.431372,
        #     "Bessell_R": 0.190278,
        #     "Kepler_Kp": 0.122666,
        #     "TESS": 0.363778,
        #     "VISTA_Z": 	0.502,
        #     "VISTA_Y": 0.600,
        #     "VISTA_J": 0.916,
        #     "VISTA_H": 1.366,
        #     "VISTA_Ks": 1.827}
        # for filt, corr in Vega_AB_Correction.items():
        #    if filt in dataframe.columns:
        #        dataframe[filt] = (dataframe[filt] + corr).round(5)
        #    else:
        #        dataframe[filt] = 99


        # for filt in Vega_AB_Correction.keys():
        #     dataframe[filt].fillna(value=99, inplace=True)

        # dataframe.fillna(value=2e-50, inplace=True)
        
        # add K-band magnitude
        k213 = dataframe["2MASS_Ks"] + 1.834505
       
        
        # reduce data frame to the needed columns
        filtlist = self.model.parms.chosen_bands
        cols = filtlist + ["mul", "mub", "Vr", "U", "V", "W", "iMass", 
                "CL", "age", "Teff", "logg", "pop", "Mass", "Mbol", "Radius", 
                "[Fe/H]","l", "b", "RA2000.0", "DEC2000.0", 
                "Dist", "x", "y", "z", "Av", "[alpha/Fe]"]
        dataframe = dataframe[cols]
        idx = dataframe.columns.get_loc("F184")+1
        dataframe.insert(idx, "K213", k213)
        
        
        for filt in filtlist:
            dataframe[filt].fillna(value=99, inplace=True)
        
        dataframe['K213'].fillna(value=99, inplace=True)
        
        dataframe.fillna(value=2e-50, inplace=True)

        if self.cat_type == "lens":
            pass

        elif self.cat_type == "source":
            # dataframe = dataframe.loc[dataframe["W146"] < 25]
            pass # upper limit should be handled in config file
            
        elif self.cat_type == "extra-bright":
            # dataframe = dataframe.loc[dataframe["W146"] < 12] #aim for ~5k+
            pass # ditto
            
        elif self.cat_type == "bright":
            # dataframe = dataframe.loc[dataframe["W146"] < 17] #aim for ~10k+
            dataframe = dataframe.loc[dataframe["W146"] > 12]

        elif self.cat_type == "mid1":
            # dataframe = dataframe.loc[dataframe["W146"] < 20]
            dataframe = dataframe.loc[dataframe["W146"] > 17] #aim for ~20k+

        elif self.cat_type == "mid2":
            # dataframe = dataframe.loc[dataframe["W146"] < 24] #aim for ~20k+
            dataframe = dataframe.loc[dataframe["W146"] > 20]

        elif self.cat_type == "faint":
            # dataframe = dataframe.loc[dataframe["W146"] < 35] #aim for ~10k+
            dataframe = dataframe.loc[dataframe["W146"] > 24]

        return dataframe
