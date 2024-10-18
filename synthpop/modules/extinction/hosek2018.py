"""
This file provides an implementation of the extinction Law from Hosek et al 2018
THIS IS A PRELIMINARY VERSION OF THE EXTINCTION LAW
"""

__all__ = ["Hosek2018", ]
__author__ = "M.J. Huston"
__date__ = "2022-07-10"
__license__ = "GPLv3"
__version__ = "1.0.0"

from ._extinction import ExtinctionLaw
from scipy import interpolate
import numpy as np


class Hosek2018(ExtinctionLaw):
    """
    Extinction law from Hosek et al 2018 :
    """

    def __init__(self, **kwargs, ):
        self.extinction_law_name = 'Hosek2018'
        self.law_ref_wavelength = 2.14

    def Alambda_Aref(self, eff_wavelength: float, R_V: float = 3.1) -> float:
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
            
        # Check wavelegth range
        #if ((min(wavelength) < 0.8) | (max(wavelength) > 2.2)):
        #    msg = 'Extinction law not defined at wavelength. Please select value between 0.8 - 2.2 microns'
        #    raise ValueError(msg)
    
        # Define extinction law (A_lambda / A_Ks) according to Hosek+17
        wave_law = np.array([0.8059, 0.962, 1.25, 1.53, 2.14, 3.545])
        A_AKs_law = np.array([9.66, 6.29, 3.56, 2.33, 1.0, 0.50])
    
        # Interpolate over the curve with cubic spline interpolation
        spline_interp = interpolate.splrep(wave_law, A_AKs_law, k=3, s=0)

        # Evaluate law at input wavelength(s)
        A_AKs_at_wave = interpolate.splev(eff_wavelength, spline_interp)
        
        # Multiply A_lambda / A_Ks * A_Ks to produce total extinction
        Alambda_AKs = A_AKs_at_wave
    
        return Alambda_AKs
    
