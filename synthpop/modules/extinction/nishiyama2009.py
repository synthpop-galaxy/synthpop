"""
Extinction law from Nishiyama et al (2009), for the direction of the Galactic
center. We apply a a spline interpolation over the measurements provided in the paper.

Valid from 1.2 to 8.0 microns.

Source DOI: 10.1088/0004-637X/696/2/1407
"""

__all__ = ["Nishiyama2009"]
__author__ = "J. Klüter, M.J. Huston"
__date__ = "2022-11-05"

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
try:
    from ._extinction import ExtinctionLaw
except ImportError:
    from _extinction import ExtinctionLaw

class Nishiyama2009(ExtinctionLaw):
    """
    Extinction law from Nishiyama et al. 2009
    Using a spline interpolation to estimate the extinction
    between 2.214 microns  and 7.60 microns
    Use Alambda prop lambda^-2.0 for lambda l below 2 microns 
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.extinction_law_name = 'Nishiyama2009'
        self.table1, self.spline = self.determine_parameters()
        self.law_ref_wavelength = 2.14
        self.min_wavelength = 1.2
        self.max_wavelength = 8.0

    @staticmethod
    def determine_parameters():
        """
         creates an interpolation to the table data.
        """
        # load Table as numpy array 
        table1 = np.array([[ 1.250e+00,     np.nan,     np.nan,  3.020e+00,  4.000e-02],
                           [ 1.630e+00,     np.nan,     np.nan,  1.730e+00,  3.000e-02],
                           [ 2.140e+00,     np.nan,     np.nan,  1.000e+00,  0.000e+00],
                           [ 3.545e+00,  2.010e+00,  4.000e-02,  5.000e-01,  1.000e-02],
                           [ 4.442e+00,  1.640e+00,  2.000e-02,  3.900e-01,  1.000e-02],
                           [ 5.675e+00,  1.560e+00,  3.000e-02,  3.600e-01,  1.000e-02],
                           [ 7.760e+00,  1.740e+00,  4.000e-02,  4.300e-01,  1.000e-02],
                           [ 1.240e+00, -5.280e-01,  1.500e-02,  2.890e+00,  8.000e-02],
                           [ 1.664e+00, -1.610e+00,  4.000e-02,  1.620e+00,  4.000e-02],
                           [ 2.164e+00,     np.nan,     np.nan,  1.012e+00,  2.000e-03]])

        # spline interpolation needs a sorted array
        table2 = table1[np.argsort(-table1[:7, 0])]
        spline = UnivariateSpline(-np.log10(table2[:, 0]), np.log10(table2[:, 3]),
            w=table2[:, 3] * 200 / np.maximum(table2[:, 4], 1e-3))

        return table1, spline

    def Alambda_Aref(self, eff_wavelength: float or np.ndarray) -> float:
        """
        Given an effective wavelength lambda_eff, calculate the extinction ratio A_lambda/A_ref

        Parameters
        ----------
        eff_wavelength : float
            wavelength to compute extinction ratio at [microns]
        """
        if eff_wavelength < 2.140:
            return (eff_wavelength/2.140) ** -2
        else:
            return 10 ** self.spline(-np.log10(eff_wavelength))
