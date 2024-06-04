"""
This file provides an implementation of the extinction Law from Damineli et al 2016.
"""

__all__ = ["Damineli", ]
__author__ = "M.J. Huston"
__date__ = "2024-04-18"
__license__ = "GPLv3"
__version__ = "1.0.0"

from ._extinction import ExtinctionLaw
import numpy as np


class Damineli(ExtinctionLaw):
    """
    Good for 0.4-4.8 microns.
    """

    def __init__(self, **kwargs, ):
        self.extinction_law_name = 'Damineli'
        self.law_ref_wavelength = 2.159

    def Alambda_AV(self, eff_wavelength: float, R_V: float = 3.1) -> float:
        return

    def Alambda_AKs(self, eff_wavelength: float) -> float:
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

        x = np.log10(2.159/eff_wavelength)
        Alam_AKs = 10**(-0.015 + 2.330*x + 0.522*x**2 - 3.001*x**3 + 2.034*x**4)

        return Alam_AKs


    def Alambda_Afilt(self, eff_wavelength: float, R_V: float = 3.1) -> float:
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
        Afilt_AKs = self.Alambda_AKs(self.ref_wavelength)
        Alam_AKs = self.Alambda_AKs(eff_wavelength)
        # Return  A_lambda   A_lambda    A_V
        #        -------- = -------- * ------
        #         A_filt      A_V      A_filt
        return Alam_AKs / Afilt_AKs




