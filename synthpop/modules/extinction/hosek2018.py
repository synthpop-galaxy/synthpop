"""
Extinction law from Hosek et al (2018), representative of highly reddened
regions like young star clusters and the Galactic center.

Valid from 0.8 to 2.2 microns.

Source DOI: 10.3847/1538-4357/aaabbb

Code also available at https://github.com/mwhosek/extlaw_H18
"""

__all__ = ["Hosek2018", ]
__author__ = "M.J. Huston"
__date__ = "2022-07-10"

try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw
from scipy import interpolate
import numpy as np


class Hosek2018(ExtinctionLaw):
    """
    Extinction law from Hosek et al. (2018)
    """

    def __init__(self, **kwargs):
        self.extinction_law_name = 'Hosek2018'
        self.law_ref_wavelength = 2.14
        self.min_wavelength = 0.8
        self.max_wavelength = 2.2

    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the relative extinction A_lambda/A_V

        Parameters
        ----------
        eff_wavelength : float
            Effective Wavelength of the filter for which the extinction should be determined.
            in micrometer
        """
    
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
    
