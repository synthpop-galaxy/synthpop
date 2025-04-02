"""
Extinction law that pulls from Surot et al. (2020) to modify the IR end of
the extinction law of O'Donnell (1994) / Cardelli et al. (1989).

Valid from: 0.25 to 3.5 microns

Reference DOIs: 10.1086/167900, 10.1086/173713, 10.1051/0004-6361/202038346

Derived by Marz Newman
"""

__all__ = ['SODC']
__author__ = "M. Newman, M.J. Huston, J. Klüter"
__date__ = "2024-05-01"

try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw
    
class SODC(ExtinctionLaw):
    """
    New combination extinction law
       
    Attributes
    -------
    R_V=3.1 : float
        optical total-to-selective extinction ratio
    alpha=2.255 : float
        power law index
    """
    
    def __init__(self, R_V: float = 3.1, alpha: float =2.255,  **kwargs):
        self.extinction_law_name = 'NewExtLaw'
        self.law_ref_wavelength = 0.549
        self.R_V = R_V
        self.alpha=alpha
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
            a = 1 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 - 1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8
            b = 1.952*y + 2.908*y**2 - 3.989*y**3 - 7.985*y**4 + 11.102*y**5 + 5.491*y**6 - 10.805*y**7 + 3.347*y**8

        else:
            a = 0.53974 * x ** self.alpha
            b = -0.495567 * x ** self.alpha
        
        #return Al_AV
        return a + b / self.R_V
