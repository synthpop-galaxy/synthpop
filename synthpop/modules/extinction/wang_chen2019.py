"""
Extinction law from WANG S. and CHEN X. (2019).
The optical to mid-infrared extinction law based on the
APOGEE, Gaia DR2, Pan-STARRS1, SDSS, APASS, 2MASS, and WISE surveys.

Valid from 0.3 to 3.3 microns

Source: 2019ApJ...877..116W - Astrophys. J., 877, 116-116 (2019/June-1)
"""

__all__ = ['WangChen2019']
__author__ = "J. Klüter, M.J. Huston"
__date__ = "2022-11-05"

try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw
    
class WangChen2019(ExtinctionLaw):
    """
    Extinction law from Wang & Chen (2019)
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.extinction_law_name = 'WangChen2019'
        self.law_ref_wavelength = 0.549
        self.min_wavelength = 0.3
        self.max_wavelength = 3.33

    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the extinction ratio A_lambda/A_ref

        Parameters
        ----------
        eff_wavelength : float
            wavelength to compute extinction ratio at [microns]
        """

        x = 1. / eff_wavelength
        if x >= 1.0:
            y = x - 1.82
            Al_Aref = (1 + 0.7499 * y - 0.1086 * y ** 2 - 0.08909 * y ** 3 + 0.02905 * y ** 4
                     + 0.01069 * y ** 5 + 0.001707 * y ** 6 - 0.001002 * y ** 7)
        else:
            Al_Aref = 0.3722 * x ** 2.07

        return Al_Aref
