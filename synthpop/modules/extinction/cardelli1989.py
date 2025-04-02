""" 
Extinction Law from Cardelli, Clayton, and Mathis (1989). 

Valid from 0.125 to 3.5 microns.

R_V-dependent function, with default R_V=3.1.

Source DOI: 10.1086/167900
"""

__all__ = ["Cardelli1989"]
__author__ = "J. Klüter, M.J. Huston"
__date__ = "2022-07-10"

try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw
    
class Cardelli1989(ExtinctionLaw):
    """
    Extinction law from Cardelli et al. (1989)
       
    Attributes
    -------
    R_V=3.1 : float
        optical total-to-selective extinction ratio
    """

    def __init__(self, R_V: float = 3.1, **kwargs):
        super().__init__(**kwargs)
        self.extinction_law_name = "Cardelli1989"
        self.law_ref_wavelength = 0.549
        self.R_V=3.1
        self.min_wavelength = 0.125
        self.max_wavelength = 3.5

    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the extinction ratio A_lambda/A_ref

        Parameters
        ----------
        eff_wavelength : float
            wavelength to compute extinction ratio at [microns]
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

        return a + b / self.R_V
