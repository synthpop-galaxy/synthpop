"""
Bulge density profile from from unpublished work related to Koshimoto et al 2021
"""

__all__ = ["Koshimoto2021BulgeB", ]
__author__ = "M.J. Huston"
__date__ = "2022-02-02"

import numpy as np
import scipy.special
from ._population_density import PopulationDensity

class Koshimoto2021BulgeB(PopulationDensity):
    """
    Bessell function bulge density distributions from Koshimoto for
    (private communication, unpublished work related to Koshimoto+21)
    
    Attributes
    ----------
    x0 : float
    y0 : float
    z0 : float
    C_par : float
    C_perp : float
    bar_ang : float
    """
        
    def __init__(
            self, x0=0.849918751795326, y0=0.339420928043361, z0=0.286256780667543, n0=7.53034e9,
            C_perp=1.28032639342739, C_par=3.24013809549932, bar_angle=27, **kwargs
            ):
        super().__init__(**kwargs)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.n0 = n0
        self.C_par = C_par
        self.C_perp = C_perp
        self.density_unit = 'mass'
        self.bar_ang = bar_angle * np.pi / 180

    def density(self, r, theta, z):
        alpha = np.pi/2.0 - self.bar_ang

        # switch to cartestian coordinates
        x = r * np.cos(theta)
        y = r * np.sin(theta)

        # Rotation to account for rotation angle
        xf = x * np.cos(alpha) + y * np.sin(alpha)
        yf = -x * np.sin(alpha) + y * np.cos(alpha)
        zf = z
        
        #calculations
        rs=((np.abs(xf / self.x0) ** self.C_perp
             + np.abs(yf / self.y0) ** self.C_perp) ** (self.C_par/self.C_perp)
            + np.abs(zf / self.z0) ** self.C_par ) ** (1/self.C_par)
        rho = self.n0 * scipy.special.kn(0, rs)

        return rho
