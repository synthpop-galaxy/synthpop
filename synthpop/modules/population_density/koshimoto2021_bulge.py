"""
Bulge density profile from Koshimoto et al. (2021)
"""

__all__ = ["Koshimoto2021Bulge", ]
__author__ = "M.J. Huston"
__date__ = "2022-02-02"

import numpy as np
import scipy.special
from ._population_density import PopulationDensity

class Koshimoto2021Bulge(PopulationDensity):
    """
    Bulge density distribution options from Koshimoto+21 (+ unpublished related work,
    N. Koshimoto, private communication)
    
    Attributes
    ----------
    parameterization : str
        parameterization options, with 'E' for Exponential, 'G' for Gaussian, 
        and 'B' for Bessel function.
    x0 : float
        scale length along x' axis
    y0 : float
        scale length along y' axis
    z0 : float
        scale length along z' axis
    C_par : float
        bar shape parameter
    C_perp : float
        bar shape parameter
    bar_angle_deg : float
        angle of the bar from the GC line of sight in degrees
    R_c : float
        cutoff radius
    X_shape : boolean
        if true, apply as an X-shaped structure
    b_X : float
        the slope of the X-shaped structure
    """
        
    def __init__(
            self, parameterization, rho0, x0, y0, z0, C_perp, C_par, R_c, 
            X_shape=False, b_X=None, bar_angle_deg=27, **kwargs
            ):
        # these were the defaults we phased out:
        #    parameterization='B', x0=0.849918751795326, y0=0.339420928043361, z0=0.286256780667543, 
        #    rho0=7.53034e9, C_perp=1.28032639342739, C_par=3.24013809549932, bar_angle_deg=27, 
        super().__init__(**kwargs)
        assert (parameterization in ['E', 'G', 'B']), f"Invalid Koshimoto2021Bulge parameterization" \
            f" '{parameterization}'. Options are 'E' (exponential), 'G' (gaussian), or 'B' (bessell)."
        self.parameterization = parameterization
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.rho0 = rho0
        self.C_par = C_par
        self.C_perp = C_perp
        self.density_unit = 'mass'
        self.bar_angle_rad = bar_angle_deg * np.pi / 180
        param_functions = {'E':self.param_function_e,
                           'G':self.param_function_g,
                           'B':self.param_function_b}
        self.param_function = param_functions[parameterization]
        self.X_shape=X_shape
        self.R_c = R_c

    @staticmethod
    def param_function_e(rs):
        return np.exp(-rs)

    @staticmethod
    def param_function_g(rs):
        return np.exp(-0.5*rs**2)

    @staticmethod
    def param_function_b(rs):
        return scipy.special.kn(0, rs)

    @staticmethod
    def cutoff_function(x):
        # Eqn 14
        return np.exp(-x**2) ** (x>0)

    def rs_function(xp, yp, zp):
        # Eqn 16
        rs = ((np.abs(xp / self.x0) ** self.C_perp
             + np.abs(yp / self.y0) ** self.C_perp) ** (self.C_par/self.C_perp)
             + np.abs(zp / self.z0) ** self.C_par ) ** (1/self.C_par)
        return rs

    def density(self, r, phi_rad, z):

        # Align coordinates with the bar,
        xp = -r * np.cos(phi_rad - self.bar_ang)
        yp = r * np.sin(phi_rad - self.bar_ang)
        zp = z
        
        # Do the math
        if not self.X_shape:
            # Eqn 13
            rs = rs_function(xp, yp, zp)
            rho = self.rho0 * self.param_function(rs) * self.cutoff_function((r-self.R_c)/0.5)
        else:
            # Eqn 17
            rs1 = rs_function(xp-self.b_X*zp, yp, zp)
            rs2 = rs_function(xp+self.b_X*zp, yp, zp)
            rho = self.rho0 * \
                    (self.param_function(rs1) + self.param_function(rs2)) * \
                    self.cutoff_function((r-self.R_c)/0.5)

        return rho
