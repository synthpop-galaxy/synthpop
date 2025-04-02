"""
Kinematic module for the based on the Koshimoto et al (2021) bulge kinematic model
"""

__all__ = ['Koshimoto2021Bulge']
__author__ = "M.J. Huston"
__date__ = "2023-05-03"

import numpy as np
from .. import const
from ._kinematics import Kinematics

class Koshimoto2021Bulge(Kinematics):
    """
    Kinematic module for the bulge based on Koshimoto et al. (2021)
    """

    def __init__(
            self, v0_stream, y0_stream, C_par_r, C_perp_r, C_par_z, C_perp_z, h0_r, h0_z, sigma_i0,
            sigma_i1, omega_p, bar_angle=27, **kwargs
            ):
        super().__init__(**kwargs) # initialises self.coord_transform & self.density_class
        self.v0_stream = v0_stream  # km/s
        self.y0_stream = y0_stream  # pc
        self.C_par_r = C_par_r  # unitless
        self.C_perp_r = C_perp_r  # unitless
        self.C_par_z = C_par_z  # unitless
        self.C_perp_z = C_perp_z  # unitless
        self.h0_r = h0_r  # kpc
        self.h0_z = h0_z  # kpc
        self.sigma_i0 = sigma_i0  # km/s
        self.sigma_i1 = sigma_i1  # km/s
        self.omega_p = omega_p  # km/s/kpc
        self.bar_ang = bar_angle*np.pi/180 #radians


    def vel_disp(self, xp, yp, zp, i):
        if i < 2:
            x0, y0, z0 = self.h0_r
            C_par, C_perp = self.C_par_r, self.C_perp_r
        else:
            x0, y0, z0 = self.h0_z
            C_par, C_perp = self.C_par_z, self.C_perp_z
        r_s = (((xp / x0) ** C_perp + (yp / y0) ** C_perp) ** (C_par / C_perp) + (
                    zp / z0) ** C_par) ** (1 / C_par)
        f_E = np.exp(-r_s)
        sigma_i = self.sigma_i0[i] + self.sigma_i1[i] * f_E
        return sigma_i

    def draw_random_velocity(self, x, y, z, **kwargs):
        """
        Generate a random u,v,w velocity vector given galactic x,y,z coordinates

        Parameters
        ----------
        x : float [kpc]
        y : float [kpc]
        z : float [kpc]

        Returns
        -------
        u : float [km/s]
        v : float [km/s]
        w : float [km/s]
        """
        # Calculate galacocentric radii
        R = np.sqrt(x ** 2 + y ** 2)
        # Calculate the angle --- i think that is unnecessary
        # theta = np.arctan2(y,x)  #arctan(y/x)
        # sigma_i = vel_disp(x,y,z)
        # Rotate to be in plane of galactic bar -> xp axis aligned with major axis of bar
        alpha=self.bar_ang

        xp = x * np.cos(alpha) + y * np.sin(alpha)
        yp = -x * np.sin(alpha) + y * np.cos(alpha)
        zp = z
        # Stream velocityy
        v_x_stream = self.v0_stream * (1 - np.exp(-(R / self.y0_stream))) * (-1) ** (
                    1 - (yp > 0).astype(int))
        # Solid body velocity
        v_y_sb = self.omega_p * R * (-1) ** (1 - (xp < 0).astype(int))

        # velocity dispersions
        sigma_x = self.vel_disp(abs(xp), abs(yp), abs(zp), 0)
        sigma_y = self.vel_disp(abs(xp), abs(yp), abs(zp), 1)
        sigma_z = self.vel_disp(abs(xp), abs(yp), abs(zp), 2)

        # Draw random deviations from circular velocity
        dvx = np.random.normal(0, sigma_x)
        dvy = np.random.normal(0, sigma_y)
        dvz = np.random.normal(0, sigma_z)

        # Calculate velocities in bar frame
        vxp = v_x_stream + dvx
        vyp = v_y_sb + dvy
        vzp = dvz

        # Convert back into Galactic frame
        rot = -alpha
        vx = -vxp * np.cos(rot) + vyp * np.sin(rot)
        vy = vxp * np.sin(rot) + vyp * np.cos(rot)
        vz = vzp

        return vx, vy, vz
