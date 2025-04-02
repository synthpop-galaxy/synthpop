"""
Velocity gradient kinematics module

This module applies solid body rotation until it reaches the LSR rotation velocity,
then uses that value. It can apply velocity dispersions in three dimensions.
"""

__all__ = ['VelocityGradient']
__author__ = "M.J. Huston"
__date__ = "2023-03-22"

from typing import Tuple
from types import ModuleType
import numpy as np
from .. import const
from ._kinematics import Kinematics
from .. import default_sun

class VelocityGradient(Kinematics):
    """
    Velocity gradient kinematics module

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

    Methods
    -------
    draw_random_velocity(self, x: ndarray, y: ndarray, z: ndarray, mass: ndarray = None,
                all_x: ndarray = None, all_y: ndarray = None, all_z: ndarray = None,
                all_mass: ndarray = None
                ) : ndarray [km/s]
        returns a random velocity of a star in km/s.
    """

    def __init__(
            self, 
            sigma_u: float, sigma_v: float, sigma_w: float,
            sun: ModuleType = None,
            vel_grad: float = 60.0,
            **kwargs
            ):
        super().__init__(**kwargs) # get sun, coord_trans and density_class
        self.kinematics_func_name = 'VelocityGradient'
        self.sun = sun if sun is not None else default_sun
        self.sigma_u = sigma_u
        self.sigma_v = sigma_v
        self.sigma_w = sigma_w
        self.vel_grad = vel_grad

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

        # Calculate rotation velocity based on inner solid body rotation and outer velocity gradient
        # Account for asymmetric drift if indicated
        rotation_velocity = np.minimum(self.vel_grad * r, self.sun.v_lsr)


        # Get velocities in the Galactic plane in the star's frame
        u1 = du
        v1 = rotation_velocity + dv

        # Rotate into Sun's frame
        u = u1 * np.cos(phi_rad) + v1 * np.sin(phi_rad)
        v = -u1 * np.sin(phi_rad) + v1 * np.cos(phi_rad)
        w = dw

        return u, v, w
