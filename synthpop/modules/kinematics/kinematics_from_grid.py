"""
Kinematic module that interpolates values from a grid
"""

__all__ = ['KinematicsFromGrid']
__author__ = "M.J. Huston"
__date__ = "2024-04-17"

from typing import Tuple
from types import ModuleType
import pandas as pd
import numpy as np
from .. import const
from ._kinematics import Kinematics
from .. import default_sun
from scipy.interpolate import LinearNDInterpolator

class KinematicsFromGrid(Kinematics):
    """
    Grid interpolation kinematic module
    
    Attributes
    ----------
    kinematics_func_name : str
        name of the Kinematics Class
    moment_file : str
        name of a file containing a grid of kinematic information
        required columns are: r, z, v_phi, sigma_phi, signa_r, sigma_z
        units must be kpc for distance, and km/s for velocity
        file must be whitespace delimited and have comments marked with '#'
    """

    def __init__(
            self, moment_file=None, sun=None,
            **kwargs
            ):
        super().__init__(**kwargs) # get sun, coord_trans and density_class
        # Open the file and create interpolators for rotational velocity and velocity dispersions
        dat = pd.read_csv(const.MOMENTS_DIR + '/' + moment_file,
            delim_whitespace=True, comment='#')
        self.interpolate_v_phi = LinearNDInterpolator(list(zip(dat['r'],dat['z'])), 
            dat['v_phi'], fill_value=0.0, rescale=False)
        self.interpolate_sigma_phi = LinearNDInterpolator(list(zip(dat['r'],dat['z'])), 
            dat['sigma_phi'], fill_value=0.0, rescale=False)
        self.interpolate_sigma_r = LinearNDInterpolator(list(zip(dat['r'],dat['z'])), 
            dat['sigma_r'], fill_value=0.0, rescale=False)
        self.interpolate_sigma_z = LinearNDInterpolator(list(zip(dat['r'],dat['z'])), 
            dat['sigma_z'], fill_value=0.0, rescale=False)
        self.kinematics_func_name = 'kinematics_from_grid'
        self.sun = sun if sun is not None else default_sun

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
        absz = np.abs(z)

        sigma_r = self.interpolate_sigma_r(list(zip(r,absz)))
        sigma_phi = self.interpolate_sigma_phi(list(zip(r,absz)))
        sigma_z = self.interpolate_sigma_z(list(zip(r,absz)))
        v_rot = self.interpolate_v_phi(list(zip(r,absz)))

        # Draw random deviations from circular velocity
        dv_r = np.random.normal(0, sigma_r)
        dv_phi = np.random.normal(0, sigma_phi)
        dv_z = np.random.normal(0, sigma_z)

        # Get velocities in the Galactic plane in the star's frame
        v_r = dv_r
        v_phi = v_rot + dv_phi

        # Rotate into Sun's frame
        u = v_r * np.cos(phi_rad) + v_phi * np.sin(phi_rad)
        v = -v_r * np.sin(phi_rad) + v_phi * np.cos(phi_rad)
        w = dv_z

        return u, v, w
