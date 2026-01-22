"""
Extinction law module which allows a user to access any extinction laws available
in the SPISEA package.
"""

__all__ = ['Spisea']
__author__ = "M.J. Huston"
__date__ = "2026-01-22"

try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw
from spisea import reddening
import numpy as np
    
class Spisea(ExtinctionLaw):
    """
    ExtinctionLaw subclass which accesses SPISEA extinction laws. Performs linear interpolation between
    the points provided by the SPISEA reddening law.

    Parameters
    ----------
    red_law_str : str
        string indicating reddening law, in format used by spisea.reddening.get_red_law() 
        (see SPISEA documentation)
    """

    def __init__(self, red_law_str, **kwargs):
        self.extinction_law_name = 'Spisea:'+red_law_str
        self.spisea_red_law = reddening.get_red_law(red_law_str)
        self.min_wavelength = self.spisea_red_law.low_lim*1e-4
        self.max_wavelength = self.spisea_red_law.high_lim*1e-4
        self.law_ref_wavelength = self.spisea_red_law.wave[0]*1e-4

    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the extinction ratio A_lambda/A_ref

        Parameters
        ----------
        eff_wavelength : float
            wavelength to compute extinction ratio at [microns]
        """
        obsc = np.interp(eff_wavelength*1e4, self.spisea_red_law.wave, self.spisea_red_law.obscuration,
                left=np.nan, right=np.nan)
        
        return obsc/self.spisea_red_law.obscuration[0]
