__all__ = ['Besancon2003']
__author__ = "M.J. Huston"
__date__ = "2023-03-22"
__license__ = "GPLv3"
__version__ = "1.0.0"

import warnings
import numpy as np
from .. import const
from ._kinematics import Kinematics
from typing import Tuple
from types import ModuleType


class Besancon2003(Kinematics):
    """
    The Kinematics base class for a Population class. The appropriate subclass is
    assigned based on the kinematics_func_kwargs through the "get_subclass" factory.

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
    vel_grad : float
        rotation velocity gradient
    disp_grad : float
        velocity dispersion gradient (dln(sigma_u^2)/dR)

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
            sigma_u: float, sigma_v: float, sigma_w: float,
            vel_grad: float = 60.0,
            disp_grad: float = 0,
            do_V_ad=True, const_V_ad=None,
            **kwargs
            ):
        """ Init """
        super().__init__(**kwargs)
        self.kinematics_func_name = 'Besancon2003'
        self.sigma_u = sigma_u
        self.sigma_v = sigma_v
        self.sigma_w = sigma_w
        self.vel_grad = vel_grad
        self.disp_grad = disp_grad
        self.do_V_ad = do_V_ad
        self.const_V_ad = const_V_ad

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

        # Apply velocity dispersion gradient
        sigma_u_r = self.sigma_u * np.exp((r - self.sun.r) * self.disp_grad / 2)
        # Draw random deviations from circular velocity
        du = np.random.normal(0, sigma_u_r)
        dv, dw = np.random.normal(0, [self.sigma_v, self.sigma_w], (*r.shape, 2)).T

        # Calculate asymmetric drift
        # For a single value, adopt the maximum of the given value
        # and the rotational velocity at the position
        if self.const_V_ad is not None:
            V_ad = np.minimum(self.const_V_ad, self.vel_grad * r)

        # Calculate asymmetric drift from equation in Robin et al, (2003) paper erratum
        elif self.do_V_ad:
            if self.density_class.density_unit == 'number':
                avg_mass = self.density_class.average_mass
            else:
                avg_mass = 1

            rho = self.density_class.density(r, phi_rad, z) * avg_mass
            drho_dr = self.density_class.gradient(r, phi_rad, z)[0] * avg_mass
            dsigma_u_dr = self.disp_grad / 2 * sigma_u_r
            V_ad = -sigma_u_r ** 2 / (2 * self.sun.v_lsr) * \
                (r / rho * drho_dr
                 + 2 * r / sigma_u_r * dsigma_u_dr
                 + (1 - self.sigma_v ** 2 / sigma_u_r ** 2)
                 + (1 - self.sigma_w ** 2 / sigma_u_r ** 2)
                 )
            V_ad = np.maximum(V_ad, 0)

        # If selected, don't apply an asymmetric drift
        else:
            V_ad = 0.0

        # Calculate rotation velocity based on inner solid body rotation and outer velocity gradient
        # Account for asymmetric drift if indicated
        rotation_velocity = np.minimum(self.vel_grad * r, self.sun.v_lsr) - V_ad

        # Get velocities in the Galactic plane in the star's frame
        u1 = du
        v1 = rotation_velocity + dv

        # Rotate into Sun's frame
        u = u1 * np.cos(phi_rad) + v1 * np.sin(phi_rad)
        v = -u1 * np.sin(phi_rad) + v1 * np.cos(phi_rad)
        w = dw

        return u, v, w