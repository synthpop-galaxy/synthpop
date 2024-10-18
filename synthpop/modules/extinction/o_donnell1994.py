"""
This file provides an implementation of the extinction Law from ODonnell 1994.
It is bases on Cardelli 1989.
"""

__all__ = ["ODonnell1994", ]
__author__ = "M.J. Huston"
__date__ = "2022-07-10"
__license__ = "GPLv3"
__version__ = "1.0.0"

from ._extinction import ExtinctionLaw


class ODonnell1994(ExtinctionLaw):
    """
    Extinction law from O'Donnell(1994) :
    For Alambda_AV. Good below 1/lambda=3 or so. Or at least smoother than Cardelli
    DO NOT GO ABOVE THIS, as it ramps up crazy fast
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.extinction_law_name = 'ODonnell1994'
        self.law_ref_wavelength=0.549

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

        x = 1. / eff_wavelength
        if x >= 1.1:
            y = x - 1.82
            a = 1 + y * (0.104 + y * (-0.609 + y * (0.701 + y * (1.137
                + y * (-1.718 + y * (-0.827 + y * (1.647 + y * (-0.505))))))))
            b = y * (1.952 + y * (2.908 + y * (-3.989 + y * (-7.985
                + y * (11.102 + y * (5.491 + y * (-10.805 + y * 3.347)))))))

        else:
            xpow = x ** 1.61
            a = 0.574 * xpow
            b = -0.527 * xpow

        return a + b / R_V
