"""
Kinematics class based on Sofue (2015) velocity curve for the Milky Way
"""
__all__ = ['Sofue2015']
__author__ = "M.J. Huston"
__date__ = "2025-02-04"
__license__ = "GPLv3"
__version__ = "1.0.0"

from typing import Tuple
from types import ModuleType
import numpy as np
from .. import const
from ._kinematics import Kinematics
from .. import default_sun
import scipy.special

class Sofue2015(Kinematics):
    """
    Attributes
    ----------
    kinematics_func_name : str
        name of the Kinematics Class
    sigma_u : float
        velocity dispersion in x direction
    sigma_v : float
        velocity dispersion in y direction
    sigma_w : float
        velocity dispersion in z direction

    Methods
    -------
    __init__(self, Population) : None
        initialize the Kinematics class
    draw_random_velocity(self, x: ndarray, y: ndarray, z: ndarray, mass: ndarray = None,
                all_x: ndarray = None, all_y: ndarray = None, all_z: ndarray = None,
                all_mass: ndarray = None
                ) : ndarray [km/s]
        returns a random velocity of a star in km/s.
    """

    def __init__(
            self, 
            sigma_u: float = 0.0, sigma_v: float = 0.0, sigma_w: float = 0.0,
            sun: ModuleType = None,
            vel_grad: float = 60.0,
            **kwargs
            ):
        """ Init """
        super().__init__(**kwargs) # get sun, coord_trans and density_class
        self.kinematics_func_name = 'Sofue2015'
        self.sun = sun if sun is not None else default_sun
        self.sigma_u = sigma_u
        self.sigma_v = sigma_v
        self.sigma_w = sigma_w
        # Bulge parameters
        self.M_b = 0.25e11 # Msun
        self.a_b = 0.87    # kpc
        self.kappa = 7.6695
        # Disk parameters
        self.M_d = 1.12e11 # Msun
        self.a_d = 5.73    # kpc
        self.sigma_d_0 = self.M_d/(2*np.pi*self.a_d**2)
        # Halo parameters
        self.h = 10.7       # kpc
        self.p_0 = 18.2e-3  # Msun/pc^3

    def draw_random_velocity(
            self, x: np.ndarray or float, y: np.ndarray or float,
            z: np.ndarray or float, **kwargs
            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate a random velocities u,v,w  by using a velocity dispersion

        Parameters
        ----------
        x, y, z : nparray, float [kpc]
            galactocentric coordinates

        Returns
        -------
        u, v, w : ndarray [km/s]
            velocity in galactocentric x,y,and z direction.
        """

        # Convert to Galactocentric coordinates
        r, phi_rad, z = self.coord_trans.xyz_to_rphiz(x, y, z)

        # Draw random deviations from circular velocity
        du, dv, dw = np.random.normal(0, [self.sigma_u, self.sigma_v, self.sigma_w],
                                      (*r.shape, 3)).T

        # Bulge contribution (Eq. 6-8)
        M_b_r = self.M_b * scipy.special.gammainc(8, self.kappa*(r/self.a_b)**(1./4))
        v_b = np.sqrt(const.G*M_b_r*const.Msun_kg/(r*1e3*const.m_per_pc)) * 1e-3
        
        # Dist contribution (Eqs 9-11; see also Eq 12 from Freeman 1970)
        r_tilde = r/self.a_d
        v_tilde = np.sqrt(0.5*r_tilde**2 *
                          (scipy.special.i0(0.5*r_tilde)*scipy.special.i1(0.5*r_tilde) +
                          scipy.special.kn(0,0.5*r_tilde)*scipy.special.kn(0,0.5*r_tilde)))
        v_d = np.sqrt(const.G * self.M_d * const.Msun_kg / (self.a_d*1e3*const.m_per_pc)) * v_tilde * 1e-3
        
        # Halo contribution (Eq. 12-14)
        r_h = r/self.h
        M_h_r = 4*np.pi * self.p_0*self.h**3 * 1e9 * (np.log(1+r_h) - r_h/(1+r_h))
        v_h = np.sqrt(const.G*M_h_r*const.Msun_kg/(r*1e3*const.m_per_pc)) * 1e-3

        # Combined rotation velocity (Eq. 5)
        rotation_velocity = np.sqrt(v_b**2 + v_d**2 + v_h**2)

        # Get velocities in the Galactic plane in the star's frame
        u1 = du
        v1 = rotation_velocity + dv

        # Rotate into Sun's frame
        u = u1 * np.cos(phi_rad) + v1 * np.sin(phi_rad)
        v = -u1 * np.sin(phi_rad) + v1 * np.cos(phi_rad)
        w = dw

        return u, v, w
