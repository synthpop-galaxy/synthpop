"""
Extinction law from Fitzpatrick & Massa (2009), where an alternative to
a lower law model allows for good fit among many lines of sight when the 
free parameters are varied appropriately.

Function dependent on R_V and alpha parameters, 
typically alpha~2.5 and R_V=3 (our defualt), or
alpha~1.8, R_V~5.

Valid from 0.5-3 microns.

Source DOI: 10.1088/0004-637X/699/2/1209
"""

__all__ = ["Fitzpatrick2009", ]
__author__ = "M.J. Huston"
__date__ = "2024-06-01"

try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw
    
class Fitzpatrick2009(ExtinctionLaw):
    """
    Extinction law from Fitzpatrick et al. (2009)
       
    Attributes
    -------
    R_V=3.1 : float
        optical total-to-selective extinction ratio
    alpha=2.5 : float
        exponent in extinction law function
    """
    
    def __init__(self,R_V: float = 3.1, alpha: float = 2.5, **kwargs):
        self.extinction_law_name = 'Fitzpatrick2009'
        self.law_ref_wavelength = 0.549
        self.R_V=R_V
        self.alpha=alpha
        self.min_wavelength = 0.5
        self.max_wavelength = 3.0

    def Alambda_Aref(self, eff_wavelength: float) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the extinction ratio A_lambda/A_ref

        Parameters
        ----------
        eff_wavelength : float
            wavelength to compute extinction ratio at [microns]
        """

        k = (0.349 + 2.087*self.R_V) * (1.0 / (1.0 + (eff_wavelength / 0.507)**self.alpha)) - self.R_V

        return (k / self.R_V) + 1.

