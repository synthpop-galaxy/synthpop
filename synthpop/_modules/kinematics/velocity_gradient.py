""" This model provides a Kinematic class using a Velocity Gradient """
__all__ = ['VelocityGradient']
__author__ = "M.J. Houston"
__date__ = "2022-07-10"
__license__ = "GPLv3"
__version__ = "1.0.0"

from typing import Tuple
import numpy as np
from .. import const
from ._kinematics import Kinematics


class VelocityGradient(Kinematics):
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
            self, sigma_u: float, sigma_v: float, sigma_w: float, vel_grad: float = 60.0,
            **kwargs
            ):
        """ Init """
        self.kinematics_func_name = 'VelocityGradient'
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

        # Calculate galactocentric radii
        r = np.sqrt(x ** 2 + y ** 2)
        # Calculate the angle
        theta = np.arctan2(y, -x)
        # Draw random deviations from circular velocity
        du, dv, dw = np.random.normal(
            0, [self.sigma_u, self.sigma_v, self.sigma_w], (*r.shape, 3)
            ).T

        # Make an array to fill in the values of circular velocity
        # Set the circular velocity based whether the gradiant is less than the radius
        # self.vel_grad * r if r < const.V_LSR/self.vel_grad else const.V_LSR
        rotation_velocity = np.minimum(self.vel_grad * r, const.V_LSR)

        # assume the star is along the line of sight between sun and galaxy
        u1 = du
        v1 = rotation_velocity + dv

        # rotate the velocity
        u = u1 * np.cos(theta) + v1 * np.sin(theta)
        v = -u1 * np.sin(theta) + v1 * np.cos(theta)
        w = dw

        return u, v, w
