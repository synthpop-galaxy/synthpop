"""
Kinematic module from the Besancon Model, as described by Robin et al. (2003), with
the rotation curve as tabulated in the Galaxia code data files (Sharma et al., 2011).

The module can be adapted for all of the populations, based
on assigning the appropriate keyword arguments for velocity dispersions, velocity dispersion
gradients, and asymmetric drift handling.
"""

__all__ = ['Besancon2003']
__author__ = "M.J. Huston"
__date__ = "2025-03-27"

import warnings
import numpy as np
from .. import const
from ._kinematics import Kinematics
from typing import Tuple
from types import ModuleType
from scipy.interpolate import interp1d

class Besancon2003(Kinematics):
    """
    Kinematic module for the Robin et al. (2003) Besancon Model.

    Attributes
    ----------
    sigma_u : float
        velocity dispersion in x direction
    sigma_v : float
        velocity dispersion in y direction
    sigma_w : float
        velocity dispersion in z direction
    disp_grad : float
        velocity dispersion gradient (dln(sigma_u^2)/dR)
    do_V_ad : boolean
        if True, calculate assymmetric drift via the equation in Robin+2003
    const_V_ad : float
        if provided, this will be used as a constant asymmetric drift value

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
        self.disp_grad = disp_grad
        self.do_V_ad = do_V_ad
        self.const_V_ad = const_V_ad
        self.grid_r = np.array([1.000e-01, 2.000e-01, 3.000e-01, 4.000e-01, 5.000e-01, 6.000e-01,
           7.000e-01, 8.000e-01, 9.000e-01, 1.000e+00, 1.100e+00, 1.200e+00,
           1.300e+00, 1.400e+00, 1.500e+00, 1.600e+00, 1.700e+00, 1.800e+00,
           1.900e+00, 2.000e+00, 2.100e+00, 2.200e+00, 2.300e+00, 2.400e+00,
           2.500e+00, 2.600e+00, 2.700e+00, 2.800e+00, 2.900e+00, 3.000e+00,
           3.100e+00, 3.200e+00, 3.300e+00, 3.400e+00, 3.500e+00, 3.600e+00,
           3.700e+00, 3.800e+00, 3.900e+00, 4.000e+00, 4.100e+00, 4.200e+00,
           4.300e+00, 4.400e+00, 4.500e+00, 4.600e+00, 4.700e+00, 4.800e+00,
           4.900e+00, 5.000e+00, 5.100e+00, 5.200e+00, 5.300e+00, 5.400e+00,
           5.500e+00, 5.600e+00, 5.700e+00, 5.800e+00, 5.900e+00, 6.000e+00,
           6.100e+00, 6.200e+00, 6.300e+00, 6.400e+00, 6.500e+00, 6.600e+00,
           6.700e+00, 6.800e+00, 6.900e+00, 7.000e+00, 7.100e+00, 7.200e+00,
           7.300e+00, 7.400e+00, 7.500e+00, 7.600e+00, 7.700e+00, 7.800e+00,
           7.900e+00, 8.000e+00, 8.100e+00, 8.200e+00, 8.300e+00, 8.400e+00,
           8.500e+00, 8.600e+00, 8.700e+00, 8.800e+00, 8.900e+00, 9.000e+00,
           9.100e+00, 9.200e+00, 9.300e+00, 9.400e+00, 9.500e+00, 9.600e+00,
           9.700e+00, 9.800e+00, 9.900e+00, 1.000e+01, 1.040e+01, 1.080e+01,
           1.120e+01, 1.160e+01, 1.200e+01, 1.240e+01, 1.280e+01, 1.320e+01,
           1.360e+01, 1.400e+01, 1.440e+01, 1.480e+01, 1.520e+01, 1.560e+01,
           1.600e+01, 1.640e+01, 1.680e+01, 1.720e+01, 1.760e+01, 1.800e+01,
           1.840e+01, 1.880e+01, 1.920e+01, 1.960e+01, 2.000e+01, 2.040e+01,
           2.080e+01, 2.120e+01, 2.160e+01, 2.200e+01, 2.240e+01, 2.280e+01,
           2.320e+01, 2.360e+01, 2.400e+01, 2.440e+01, 2.480e+01, 2.520e+01,
           2.560e+01, 2.600e+01, 2.640e+01, 2.680e+01, 2.720e+01, 2.760e+01,
           2.800e+01, 2.840e+01, 2.880e+01, 2.920e+01, 2.960e+01, 3.000e+01,
           3.040e+01, 3.080e+01, 3.120e+01, 3.160e+01, 3.200e+01, 3.240e+01,
           3.280e+01, 3.320e+01, 3.360e+01, 3.400e+01, 3.440e+01, 3.480e+01,
           3.520e+01, 3.560e+01, 3.600e+01, 3.640e+01, 3.680e+01, 3.720e+01,
           3.760e+01, 3.800e+01, 3.840e+01, 3.880e+01, 3.920e+01, 3.960e+01,
           4.000e+01, 4.040e+01, 4.080e+01, 4.120e+01, 4.160e+01, 4.200e+01,
           4.240e+01, 4.280e+01, 4.320e+01, 4.360e+01, 4.400e+01, 4.440e+01,
           4.480e+01, 4.520e+01, 4.560e+01, 4.600e+01, 4.640e+01, 4.680e+01,
           4.720e+01, 4.760e+01, 4.800e+01, 4.840e+01, 4.880e+01, 4.920e+01,
           4.960e+01, 5.000e+01, 5.450e+01, 5.900e+01, 6.350e+01, 6.800e+01,
           7.250e+01, 7.700e+01, 8.150e+01, 8.600e+01, 9.050e+01, 9.500e+01,
           9.950e+01, 1.040e+02, 1.085e+02, 1.130e+02, 1.175e+02, 1.220e+02,
           1.265e+02, 1.310e+02, 1.355e+02, 1.400e+02, 1.445e+02, 1.490e+02,
           1.535e+02, 1.580e+02, 1.625e+02, 1.670e+02, 1.715e+02, 1.760e+02,
           1.805e+02, 1.850e+02, 1.895e+02, 1.940e+02, 1.985e+02, 2.030e+02,
           2.075e+02, 2.120e+02, 2.165e+02, 2.210e+02, 2.255e+02, 2.300e+02,
           2.345e+02, 2.390e+02, 2.435e+02, 2.480e+02, 2.525e+02, 2.570e+02,
           2.615e+02, 2.660e+02, 2.705e+02, 2.750e+02, 2.795e+02, 2.840e+02,
           2.885e+02, 2.930e+02, 2.975e+02, 3.020e+02, 3.065e+02, 3.110e+02,
           3.155e+02, 3.200e+02, 3.245e+02, 3.290e+02, 3.335e+02, 3.380e+02,
           3.425e+02, 3.470e+02, 3.515e+02, 3.560e+02, 3.605e+02, 3.650e+02,
           3.695e+02, 3.740e+02, 3.785e+02, 3.830e+02, 3.875e+02, 3.920e+02,
           3.965e+02, 4.010e+02, 4.055e+02, 4.100e+02, 4.145e+02, 4.190e+02,
           4.235e+02, 4.280e+02, 4.325e+02, 4.370e+02, 4.415e+02, 4.460e+02,
           4.505e+02, 4.550e+02, 4.595e+02, 4.640e+02, 4.685e+02, 4.730e+02,
           4.775e+02, 4.820e+02, 4.865e+02, 4.910e+02, 4.955e+02, 5.000e+02])
        self.grid_vcirc = np.array([742.0076 , 524.89887, 429.04287, 372.29891, 334.04801, 306.30727,
           285.28507, 268.8699 , 255.82296, 245.30855, 236.78942, 229.84636,
           224.20118, 219.60835, 215.90561, 212.93254, 210.58738, 208.75792,
           207.37662, 206.36207, 205.66848, 205.23457, 205.02935, 205.0063 ,
           205.14398, 205.40687, 205.78074, 206.23774, 206.76874, 207.35228,
           207.98243, 208.6424 , 209.32913, 210.02914, 210.74104, 211.45449,
           212.16915, 212.87666, 213.57787, 214.26603, 214.9423 , 215.60155,
           216.24531, 216.86923, 217.47529, 218.06005, 218.62531, 219.16847,
           219.69151, 220.1921 , 220.67227, 221.13034, 221.56804, 221.98403,
           222.38024, 222.75538, 223.11124, 223.44703, 223.76426, 224.06223,
           224.34261, 224.60475, 224.85006, 225.07824, 225.29058, 225.48668,
           225.66794, 225.83408, 225.98617, 226.12418, 226.24918, 226.36096,
           226.4606 , 226.54811, 226.62427, 226.68915, 226.74365, 226.78769,
           226.82207, 226.84695, 226.86293, 226.87007, 226.86916, 226.86019,
           226.84373, 226.82   , 226.78949, 226.75225, 226.7089 , 226.65953,
           226.60452, 226.54408, 226.47865, 226.40823, 226.33333, 226.25406,
           226.1707 , 226.08341, 225.99258, 225.89909, 225.49482, 225.04726,
           224.57008, 224.07196, 223.56024, 223.04105, 222.51948, 221.99971,
           221.48521, 220.97874, 220.48251, 219.99825, 219.52731, 219.07066,
           218.62899, 218.20274, 217.79218, 217.39735, 217.01818, 216.65449,
           216.30602, 215.97238, 215.65319, 215.34799, 215.05629, 214.77758,
           214.51136, 214.25711, 214.01429, 213.78239, 213.56092, 213.34936,
           213.14724, 212.95409, 212.76947, 212.59294, 212.42408, 212.26251,
           212.10785, 211.95974, 211.81783, 211.68182, 211.55139, 211.42625,
           211.30613, 211.19078, 211.07995, 210.97339, 210.87091, 210.7723 ,
           210.67735, 210.58589, 210.49776, 210.41278, 210.33079, 210.25166,
           210.17526, 210.10146, 210.03011, 209.96113, 209.8944 , 209.82981,
           209.76727, 209.70669, 209.64799, 209.59107, 209.53587, 209.48231,
           209.43032, 209.37983, 209.33079, 209.28314, 209.23681, 209.19175,
           209.14792, 209.10527, 209.06375, 209.02331, 208.98392, 208.94554,
           208.90813, 208.87165, 208.83607, 208.80136, 208.76748, 208.73441,
           208.70212, 208.67057, 208.63976, 208.60965, 208.58021, 208.55142,
           208.52327, 208.49574, 208.4688 , 208.44243, 208.41662, 208.39134,
           208.36659, 208.34357, 208.11273, 207.91287, 207.74732, 207.60787,
           207.48873, 207.38572, 207.29576, 207.21648, 207.14608, 207.08313,
           207.0265 , 206.97528, 206.92872, 206.8862 , 206.84722, 206.81135,
           206.77823, 206.74755, 206.71905, 206.6925 , 206.6677 , 206.64449,
           206.62271, 206.60224, 206.58295, 206.56475, 206.54754, 206.53125,
           206.5158 , 206.50112, 206.48716, 206.47387, 206.46119, 206.44909,
           206.43752, 206.42645, 206.41584, 206.40567, 206.39591, 206.38653,
           206.37751, 206.36883, 206.36047, 206.3524 , 206.34462, 206.33711,
           206.32985, 206.32283, 206.31604, 206.30946, 206.30309, 206.29691,
           206.29092, 206.2851 , 206.27945, 206.27396, 206.26861, 206.26342,
           206.25836, 206.25344, 206.24864, 206.24396, 206.2394 , 206.23494,
           206.23059, 206.22635, 206.22219, 206.21814, 206.21417, 206.2103 ,
           206.20651, 206.20279, 206.19915, 206.19557, 206.19207, 206.18863,
           206.18527, 206.18196, 206.17873, 206.17555, 206.17242, 206.16935,
           206.16632, 206.16334, 206.16041, 206.15752, 206.15469, 206.1519 ,
           206.14916, 206.14646, 206.14379, 206.14115, 206.13854, 206.13596,
           206.13341, 206.1309 , 206.12843, 206.126  , 206.12359, 206.12117])
        self.interp_vcirc = interp1d(self.grid_r, self.grid_vcirc,
                                     bounds_error=False, fill_value=(self.grid_vcirc[0], self.grid_vcirc[-1]))

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
        
        # Get velocity in rotation curve
        rotation_velocity = self.interp_vcirc(r)

        # Apply velocity dispersion gradient
        sigma_u_r = self.sigma_u * np.exp((r - self.sun.r) * self.disp_grad / 2)
        # Draw random deviations from circular velocity
        du = np.random.normal(0, sigma_u_r)
        dv, dw = np.random.normal(0, [self.sigma_v, self.sigma_w], (*r.shape, 2)).T

        # Calculate asymmetric drift
        # For a single value, adopt the maximum of the given value
        # and the rotational velocity at the position
        if self.const_V_ad is not None:
            V_ad = np.minimum(self.const_V_ad, rotation_velocity)

        # Calculate asymmetric drift from equation in Robin et al, (2003) paper erratum
        elif self.do_V_ad:
            if self.density_class.density_unit == 'number':
                avg_mass = self.density_class.average_mass
            else:
                avg_mass = 1

            rho = self.density_class.density(r, phi_rad, z) * avg_mass
            drho_dr = self.density_class.gradient(r, phi_rad, z)[0] * avg_mass
            dsigma_u_dr = self.disp_grad / 2 * sigma_u_r
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="divide by zero encountered in multiply")
                warnings.filterwarnings("ignore", message="invalid value encountered in multiply")
                rho_component = np.nan_to_num(r / rho * drho_dr, nan=0.0, posinf=0.0, meginf=0.0)
            V_ad = -sigma_u_r ** 2 / (2 * self.sun.v_lsr) * \
                (rho_component
                 + 2 * r / sigma_u_r * dsigma_u_dr
                 + (1 - self.sigma_v ** 2 / sigma_u_r ** 2)
                 + (1 - self.sigma_w ** 2 / sigma_u_r ** 2)
                 )
            V_ad = np.minimum(np.maximum(V_ad, 0), rotation_velocity)

        # If selected, don't apply an asymmetric drift
        else:
            V_ad = 0.0
            
        # Account for asymmetric drift if indicated
        rotation_velocity -= V_ad

        # Get velocities in the Galactic plane in the star's frame
        u1 = du
        v1 = rotation_velocity + dv

        # Rotate into Sun's frame
        u = u1 * np.cos(phi_rad) + v1 * np.sin(phi_rad)
        v = -u1 * np.sin(phi_rad) + v1 * np.cos(phi_rad)
        w = dw

        return u, v, w

    def mean_galactic_uvw(
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
        
        # Get velocity in rotation curve
        rotation_velocity = self.interp_vcirc(r)

        # Apply velocity dispersion gradient
        sigma_u_r = self.sigma_u * np.exp((r - self.sun.r) * self.disp_grad / 2)
        # Draw random deviations from circular velocity
        du = np.random.normal(0, sigma_u_r)
        dv, dw = np.random.normal(0, [self.sigma_v, self.sigma_w], (*r.shape, 2)).T

        # Calculate asymmetric drift
        # For a single value, adopt the maximum of the given value
        # and the rotational velocity at the position
        if self.const_V_ad is not None:
            V_ad = np.minimum(self.const_V_ad, rotation_velocity)

        # Calculate asymmetric drift from equation in Robin et al, (2003) paper erratum
        elif self.do_V_ad:
            if self.density_class.density_unit == 'number':
                avg_mass = self.density_class.average_mass
            else:
                avg_mass = 1

            rho = self.density_class.density(r, phi_rad, z) * avg_mass
            drho_dr = self.density_class.gradient(r, phi_rad, z)[0] * avg_mass
            dsigma_u_dr = self.disp_grad / 2 * sigma_u_r
            rho_component = np.nan_to_num(r / rho * drho_dr, nan=0.0)
            V_ad = -sigma_u_r ** 2 / (2 * self.sun.v_lsr) * \
                (rho_component
                 + 2 * r / sigma_u_r * dsigma_u_dr
                 + (1 - self.sigma_v ** 2 / sigma_u_r ** 2)
                 + (1 - self.sigma_w ** 2 / sigma_u_r ** 2)
                 )
            V_ad = np.minimum(np.maximum(V_ad, 0), rotation_velocity)

        # If selected, don't apply an asymmetric drift
        else:
            V_ad = 0.0
            
        # Account for asymmetric drift if indicated
        rotation_velocity -= V_ad

        # Get velocities in the Galactic plane in the star's frame
        u1 = 0
        v1 = rotation_velocity

        # DON'T rotate into Sun's frame
        u = rotation_velocity*0.0
        v = rotation_velocity
        w = rotation_velocity*0.0

        return u, v, w
