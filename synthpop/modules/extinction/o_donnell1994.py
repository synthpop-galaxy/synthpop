"""
Extinction Law from O'Donnell (1994), which builds from the
Cardelli et al. (1989) law. 

Valid from 0.25 to 3.5 microns

Source DOI: 10.1086/173713
"""

__all__ = ["ODonnell1994", ]
__author__ = "M.J. Huston"
__date__ = "2022-07-10"

try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw

class ODonnell1994(ExtinctionLaw):
    """
    Extinction law from O'Donnell (1994)
       
    Attributes
    -------
    R_V=3.1 : float
        optical total-to-selective extinction ratio
    """
    
    def __init__(self, R_V: float = 3.1, **kwargs):
        super().__init__(**kwargs)
        self.extinction_law_name = 'ODonnell1994'
        self.law_ref_wavelength=0.549
        self.R_V=R_V
        self.min_wavelength = 0.25
        self.max_wavelength = 3.5

    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the extinction ratio A_lambda/A_ref

        Parameters
        ----------
        eff_wavelength : float
            wavelength to compute extinction ratio at [microns]
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

        return a + b / self.R_V
