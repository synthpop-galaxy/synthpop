"""
Kinematic module for the based on the Koshimoto et al (2021) disk kinematic model
"""

__all__ = ['Koshimoto2021Disk']
__author__ = "M.J. Huston"
__date__ = "2023-05-03"

from typing import Tuple
from types import ModuleType
import numpy as np
from .. import const
from ._kinematics import Kinematics
from scipy.special import gamma, gammasgn, gammaln

class Koshimoto2021Disk(Kinematics):
    """
    Kinematic module for the disk based on Koshimoto et al. (2021)
    """

    def __init__(
            self,
            sigma_r_sun, sigma_z_sun,
            beta_r, beta_z,
            R_sigma_r, R_sigma_z,
            pop_age,
            **kwargs
            ):
        super().__init__(**kwargs)
        self.kinematics_func_name = 'Koshimoto2021Disk'
        self.sigma_r_sun = sigma_r_sun
        self.sigma_z_sun = sigma_z_sun
        self.beta_r = beta_r
        self.beta_z = beta_z
        self.R_sigma_r = R_sigma_r
        self.R_sigma_z = R_sigma_z
        self.pop_age = pop_age
        
        # Rotation curve from Bland-Hawthorn & Gerhard (2016)
        # [[galactocentric distance (kpc)],[rotational velocity (km/s)]]
        self.rot_curve = np.array([[0.0000000e+00, 9.4590000e-02, 1.6216000e-01, 3.5580000e-01,
        3.8468000e-01, 4.1356000e-01, 4.4244000e-01, 4.7854000e-01,
        5.1464000e-01, 5.5074000e-01, 5.9407000e-01, 6.3739000e-01,
        6.8071000e-01, 7.2403000e-01, 7.7457000e-01, 8.3233000e-01,
        9.2434000e-01, 1.0037700e+00, 1.0913400e+00, 1.2275900e+00,
        1.4414200e+00, 1.6608100e+00, 1.9640600e+00, 2.2600900e+00,
        2.5489000e+00, 2.8449300e+00, 3.1481800e+00, 3.4514300e+00,
        3.7546800e+00, 3.9785100e+00, 4.1756800e+00, 4.4189400e+00,
        4.6594700e+00, 5.0254400e+00, 5.3591800e+00, 5.6475400e+00,
        5.9496300e+00, 6.2528800e+00, 6.5561300e+00, 6.8593800e+00,
        7.1626300e+00, 7.4658800e+00, 7.7691300e+00, 8.0723800e+00,
        8.3756300e+00, 8.6788800e+00, 8.9821300e+00, 9.1951300e+00,
        9.6332700e+00, 9.8630000e+00, 9.8910500e+00, 1.5247700e+01,
        2.0205200e+01, 2.5011000e+01],
       [0.0000000e+00, 2.6168220e+01, 4.1121500e+01, 8.2262400e+01,
        8.7505510e+01, 9.2873460e+01, 9.8283010e+01, 1.0462191e+02,
        1.1103709e+02, 1.1721647e+02, 1.2304215e+02, 1.2882622e+02,
        1.3464001e+02, 1.4043596e+02, 1.4630325e+02, 1.5187448e+02,
        1.5776319e+02, 1.6350981e+02, 1.6911703e+02, 1.7503895e+02,
        1.8137693e+02, 1.8510346e+02, 1.9026516e+02, 1.9500092e+02,
        2.0044647e+02, 2.0425092e+02, 2.0746099e+02, 2.1074834e+02,
        2.1491548e+02, 2.1952495e+02, 2.2406542e+02, 2.2835254e+02,
        2.3150470e+02, 2.3417276e+02, 2.3578553e+02, 2.3745772e+02,
        2.3844409e+02, 2.3949033e+02, 2.4009074e+02, 2.4054847e+02,
        2.4034992e+02, 2.4019774e+02, 2.3993618e+02, 2.3934766e+02,
        2.3823009e+02, 2.3704757e+02, 2.3556278e+02, 2.3462386e+02,
        2.3266982e+02, 2.3166169e+02, 2.3149119e+02, 2.1587470e+02,
        2.1248535e+02, 2.1019041e+02]])

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
        # Set parameters from Koshimoto+21
        T_min, T_max = 0.01, 10
        R_d = 2.6
        c1,c2,c3,c4 = 3.822, 0.524, 0.00567, 2.13

        # Calculate velocity dispersions in directions of non-rotation
        sigma_r = self.sigma_r_sun * ((self.pop_age+T_min)/(T_max+T_min))**self.beta_r \
                  * np.exp(-(r-self.sun.r)/self.R_sigma_r)
        sigma_z = self.sigma_z_sun * ((self.pop_age+T_min)/(T_max+T_min))**self.beta_z \
                  * np.exp(-(r-self.sun.r)/self.R_sigma_z)
        
        # Set up R_g array
        R_g_arr = self.rot_curve[0]
        R_g_arr[0] = 1.0e-3
        v_c_arr = self.rot_curve[1]
        v_c_arr[0] = R_g_arr[0]/R_g_arr[1]*v_c_arr[1]
        # Get a array: r-direction velocity dispersion to circular velocity ratio at R_g positions
        a_arr = (self.sigma_r_sun * ((self.pop_age+T_min)/(T_max+T_min))**self.beta_r
                 * np.exp(-(R_g_arr-self.sun.r)/self.R_sigma_r))/v_c_arr
        a0 = self.sigma_r_sun/self.sun.v_lsr
        # Get s array: description after equation 9 in Koshimoto+21
        s_input_arr = R_g_arr/(c1*R_d*(1+R_d/self.R_sigma_r/c2))
        s_arr = 31.53*np.exp(-s_input_arr/0.2743)*((s_input_arr/0.6719)**2-1)
        # Get g array: description after equation 8 in Koshimoto+21
        # Use weird ln/exp stuff to prevent overflow
        g_input_arr = 1/(2*a_arr**2)
        g_arr = gammasgn(g_input_arr-0.5)*np.exp(g_input_arr + gammaln(g_input_arr - 0.5) - g_input_arr*np.log(g_input_arr) + 0.5*np.log(g_input_arr)) / 2
        # Get surface density for R_g array values
        surf_dens_arr = np.exp(-R_g_arr/R_d) - c3*a0**c4/R_d**2 * s_arr
        
        # Get R_g for each r based on probability distribution
        def draw_R_g(r_i):
            pdf_arr_0 = np.nan_to_num(4*np.pi**2*surf_dens_arr/g_arr * 
                    np.exp((2*np.log(R_g_arr/r_i) + 1 - R_g_arr**2/r_i**2)/(2*a_arr**2)))
            pdf_arr = pdf_arr_0*(1-(pdf_arr_0<0).astype(int))
            if sum(pdf_arr)==0:
                #print('too near GC', r)
                return r_i
            else:
                cdf_arr = np.cumsum(pdf_arr)/np.sum(pdf_arr)
                return np.interp(np.random.uniform(), cdf_arr, R_g_arr)
        if type(r) is float:
            R_g = draw_R_g(r)
        else:
            R_g = np.array(list(map(lambda r_i: draw_R_g(r_i), r)))
        
        # Calculate velocities in r, phi, z frame
        v_r = np.random.normal(0, sigma_r)
        v_phi = np.interp(R_g, self.rot_curve[0],self.rot_curve[1]) * R_g/r / (1 + 0.0374*np.abs(z)**1.34)
        v_z = np.random.normal(0, sigma_z)

        # Rotate into galactic frame
        #u = v_r * np.cos(phi_rad) + v_phi * np.sin(phi_rad)
        #v = -v_r * np.sin(phi_rad) + v_phi * np.cos(phi_rad)
        #w = v_z
        u = -v_r * np.cos(phi_rad) + v_phi * np.sin(phi_rad)
        v = v_r * np.sin(phi_rad) + v_phi * np.cos(phi_rad)
        w = v_z

        return u, v, w
