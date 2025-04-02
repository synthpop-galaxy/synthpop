"""
Extinction Law from Damineli et al (2016), derived to be representative
of the Galactic plane, |b|<5 deg.

Valid from 0.4-4.8 microns.

Source DOI: 10.1093/mnras/stw2122
"""

__all__ = ["Damineli2016", ]
__author__ = "M.J. Huston"
__date__ = "2024-04-18"

import numpy as np
try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw

class Damineli2016(ExtinctionLaw):
    """
    Extinction law from Damineli et al. (2016)
    """

    def __init__(self, **kwargs):
        self.extinction_law_name = 'Damineli2016'
        self.law_ref_wavelength = 2.159
        self.min_wavelength = 0.4
        self.max_wavelength = 4.8

    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the extinction ratio A_lambda/A_ref

        Parameters
        ----------
        eff_wavelength : float
            wavelength to compute extinction ratio at [microns]
        """

        x = np.log10(2.159/eff_wavelength)
        Alam_AKs = 10**(-0.015 + 2.330*x + 0.522*x**2 - 3.001*x**3 + 2.034*x**4)

        return Alam_AKs
