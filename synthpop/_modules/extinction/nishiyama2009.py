"""
Extinction law from Nishiyama 2009
It uses a spline interpolation form the measurements provided in the paper.
"""

__all__ = ["Nishiyama2009"]
__author__ = "J. Klüter"
__date__ = "2022-11-05"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from ._extinction import ExtinctionLaw, EXTINCTION_DIR


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

    @staticmethod
    def determine_parameters():
        """
         creates an interpolation to the table data.
        """
        # load Table as numpy array 
        table1 = np.loadtxt(f"{EXTINCTION_DIR}/nishiyama2009_table1.csv", dtype=float,
            usecols=[1, 2, 3, 4, 5], delimiter=',')
        # spline interpolation needs a sorted array
        table2 = table1[np.argsort(-table1[:7, 0])]
        spline = UnivariateSpline(-np.log10(table2[:, 0]), np.log10(table2[:, 3]),
            w=table2[:, 3] * 200 / np.maximum(table2[:, 4], 1e-3))

        # use a polynomial to fit the data        
        # pp = self.table1[:,0]>1.5
        # pp[-1]=False 
        # self.parm2 = np.polyfit(
        #           -np.log10(self.table1[pp,0]), np.log10(self.table1[pp,3]),
        #           4, w=1/np.maximum(self.table1[pp,4], 1e-5))

        return table1, spline

    def Alambda_AV(self, eff_wavelength: float or np.ndarray, R_V: float = 3.1) -> float:
        """
        Estimates the Extinction
        returns linear extrapolation for lambda < 2.140
        returns spline interpolation lambda > 2.214
        note that it is given relative to A_K 
        """
        if eff_wavelength < 2.140:
            return x ** -2
        else:
            return 10 ** self.spline(-np.log10(eff_wavelength))

    def plot_parameter(self):
        x2 = np.linspace(1, 10, 1000)
        x = np.linspace(2.140, 1, 1000)
        plt.errorbar(1 / self.table1[:, 0], self.table1[:, 3],
            yerr=self.table1[:, 4], fmt='.', label="data points")
        plt.loglog(1 / x, x ** -2 * 2.140 ** 2, label='lambda ** -2')
        plt.loglog(1 / x2, self.Alambda_AV(x2), label='spline interpolate')
        # plt.loglog(10**xx, 10**np.poly1d(self.parm2)(xx), label = 'parable interpolate')
        plt.legend()
        plt.ylabel(r"$A_{\lambda}/A_{Ks}$")
        plt.xlabel(r"$\lambda^{-1}\ [\mu met^{-1}]$")
        plt.show()
