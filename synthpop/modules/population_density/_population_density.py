"""
This file contains the base class for the population density distributions.
"""
__all__ = ["PopulationDensity", ]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-07-12"

from abc import ABC, abstractmethod
from types import ModuleType
from typing import Tuple
import numpy as np

from .. import const, default_sun
from ... import synthpop_utils as sp_utils
import pdb
from scipy.integrate import simpson

class PopulationDensity(ABC):
    """
    This is the Parent Class for population densities functions

    Attributes
    ----------
    population_density_name : str
        name of the population density
    density_unit : str
        specifier if density profile returns "number"-, "mass"- or "init_mass"-density
    (more attributes are specified in the subclasses)

    Methods
    -------
    __init__() :
        Initialize the PopulationDensity class
    density(r: ndarray, theta, ndarray, z: ndarray) : ndarray
        estimate the density at (r,theta,z)
        (specified in the subclasses)
    """
    # placeholder for average mass, and emass/imass correction
    # are set when generating a field
    average_mass  = None
    average_mass_coor = None

    def __init__(self,
            sun: ModuleType = None,
            coord_trans: ModuleType = None,
            gamma_flare: float = None,
            h_flare: float = None,
            radius_flare: float = 0,
            logger: ModuleType = None,
            **kwargs):
        """
        Initialize the Population Density class
        SubClasses MUST define a density_unit!

        Parameters
        ----------
        sun : SunInfo
            location and velocity of the sun and lsr
            see synthpop_utils/sun_info
        coord_trans: ModuleType
            the coordinate transformation package
            see synthpop_utils/coordinate_transformations
        gamma_flare, radius_flare: float
            parameters to implement the flare of the milky way

        """
        # sun sun sun, here it comes
        self.logger = logger
        self.sun = sun if sun is not None else default_sun

        self.population_density_name = 'None'
        self.density_unit = "one of 'mass', 'init_mass', or 'number'"
        self.coord_trans = coord_trans

        self.gamma_flare = 0 if gamma_flare is None and h_flare is None else gamma_flare
        self.h_flare = h_flare
        self.radius_flare = radius_flare

        # define coordinates of center of the field in degrees and radians
        self.l_deg = None
        self.l_rad = None
        self.b_deg = None
        self.b_rad = None
        # define the field
        self.field_shape = None
        self.lb_radius_deg = None
        self.l_length_deg = None
        self.b_length_deg = None

    @abstractmethod
    def density(self, r: np.ndarray, phi_rad: np.ndarray, z: np.ndarray) -> np.ndarray:
        """

        Estimates the density at the given position

        Parameters
        ----------
        r : ndarray ['kpc']
            Distance to z axis
        phi_rad : ndarray ['rad']
            azimuth angle of the stars. phi_rad = 0 is pointing towards sun.
        z : height above the galactic plane (corrected for warp of the galaxy)

        Returns
        -------
        rho : ndarray [M_sun/kpc^3 or #/kpc^3]
            density at the given location, either in number density evolved
            mass density or initial mass density should be specified in density_unit.

        """
        raise NotImplementedError('No density profile is implemented!')

    def get_kappa_flare(self,
            r_kpc: np.ndarray or float,
            gamma_flare: float = None,
            h_flare: float = None, radius_flare: float = None
            ) -> np.ndarray or float:
        """
        Estimates the correction factor for the Warp.
        The scale height should then be multiplied by kappa_flare

        Parameters
        ----------
        r_kpc : float or ndarray
            galactocentric radii
        gamma_flare: float or None
            slope of the flare
            if None use the default value from const

        radius_flare: float or None
            radius when the flare starts
            if None use the default value from const

        Returns
        -------
        kappa_flare : ndarray:
            correction factor for the scale height
        """

        if gamma_flare is None and h_flare is None:
            gamma_flare = self.gamma_flare
            h_flare = self.h_flare

        if radius_flare is None:
            radius_flare = self.radius_flare

        if gamma_flare is not None:
            return 1 + gamma_flare * np.maximum(r_kpc - radius_flare, 0)

        return np.exp(np.maximum(r_kpc - radius_flare, 0)/h_flare)

    def gradient(self, r_kpc: np.ndarray or float,
                 phi_rad: np.ndarray or float,
                 z_kpc: np.ndarray or float,
                 eps: Tuple[float] = (1e-3, 1e-4, 1e-3)
                 ) -> Tuple[np.ndarray or float, np.ndarray or float, np.ndarray or float]:
        """
        return the gradient at the given location

        Parameters
        ----------
        r_kpc :  float, ndarray [kpc]
            Radius in kpc
        phi_rad : float, ndarray [rad]
            polar angle follows the rotation of the galaxy
            zero point is at the sun.
        z_kpc : float, ndarray [kpc]
            height above/below the galactic plane
        eps: Tuple(float,float,float)
            difference used to estimate the gradient

        Returns
        -------
        dRho_dR :  float, ndarray
            Gradient in Radius
        dRho_dPhi : float, ndarray
            Gradient in polar angle direction
        dRho_dz : float, ndarray
            Gradient in z direction
        """

        dRho_dR = (self.density(r_kpc + eps[0], phi_rad, z_kpc)
                   - self.density(r_kpc - eps[0], phi_rad, z_kpc)) / (2 * eps[0])
        dRho_dPhi = (self.density(r_kpc, phi_rad + eps[1], z_kpc)
                     - self.density(r_kpc, phi_rad - eps[1], z_kpc)) / (2 * eps[1])
        dRho_dz = (self.density(r_kpc, phi_rad, z_kpc + eps[2])
                   - self.density(r_kpc, phi_rad, z_kpc - eps[2])) / (2 * eps[2])

        return dRho_dR, dRho_dPhi, dRho_dz

    def update_location(self, l_deg: float, b_deg: float, field_shape: str, field_scale_deg: float,
                        max_distance: float):
        """
        Set the location and solid_angle

        Parameters
        ----------
        l_deg, b_deg : float [deg]
            galactic longitude and latitude in degrees
        solid_angle : float [sr]
            size of the cone
        """
        self.l_deg = l_deg
        self.l_rad = l_deg * np.pi / 180.
        self.b_deg = b_deg
        self.b_rad = b_deg * np.pi / 180.
        self.field_shape = field_shape
        #  convert solid angle to half cone angle using wiki formula:
        if self.field_shape=='circle':
            self.lb_radius_deg = field_scale_deg
        elif self.field_shape=='box':
            if hasattr(field_scale_deg, 'len'):
                self.l_length_deg = field_scale_deg[0]
                self.b_length_deg = field_scale_deg[1]
            else:
                self.l_length_deg = field_scale_deg
                self.b_length_deg = field_scale_deg
        self.max_distance = max_distance

        self.logger.info("setting up density grid")
        if self.field_shape == 'circle':
            self.density_grid_d_pts = np.linspace(0, self.max_distance, 1001)
            self.density_grid_st_dir = np.linspace(0, 2 * np.pi, 101)
            self.density_grid_st_rad = np.flip(np.arccos(np.linspace(np.cos(self.lb_radius_deg*np.pi/180),1, 26)))
            d_grid, st_dir_grid, st_rad_grid = np.meshgrid(self.density_grid_d_pts, self.density_grid_st_dir, 
                                                 self.density_grid_st_rad)
            grid_shape = d_grid.shape

            # Get into physically meaningful coordinates
            delta_l_rad = st_rad_grid * np.sin(st_dir_grid)
            delta_b_rad = st_rad_grid * np.cos(st_dir_grid)
            l_rad, b_rad = self.rotate_00_to_lb(delta_l_rad.ravel(), delta_b_rad.ravel())
            l_grid = l_rad * 180 / np.pi
            if np.abs(self.l_deg)<90:
                l_grid -= 360*(l_grid>180)
            b_grid = b_rad * 180 / np.pi
            r_flat, phi_flat, z_flat = self.coord_trans.dlb_to_rphiz(d_grid.ravel(), 
                                                l_grid.ravel(), b_grid.ravel())
            # Get density at points and integrate
            self.density_grid = self.density(r_flat, phi_flat, z_flat).reshape(grid_shape)
            vol_elem = d_grid**2 * np.sin(st_rad_grid)
            self.density_int_st_rad = simpson(self.density_grid*vol_elem, x=self.density_grid_st_rad, axis=2)
            self.density_int_st_dir = simpson(self.density_int_st_rad, x=self.density_grid_st_dir, axis=0)
            self.total_mass = simpson(self.density_int_st_dir, x=self.density_grid_d_pts, axis=0)

        elif self.field_shape == 'box':
            # Box should be easier, right ?? RIGHT ???????
            # Let's make a grid
            self.density_grid_d_pts = np.linspace(0, self.max_distance, 1001)
            delta_l_rad = self.l_length_deg/2 * np.pi/180 * np.linspace(-1, 1, 51)
            delta_b_rad = self.b_length_deg/2 * np.pi/180 * np.linspace(-1, 1, 51)
            l_rad_pts, b_rad_pts = self.rotate_00_to_lb(delta_l_rad, delta_b_rad)
            self.density_grid_l_pts = l_rad_pts * 180/np.pi
            if np.abs(self.l_deg)<90:
                self.density_grid_l_pts -= 360*(self.density_grid_l_pts>180)
            self.density_grid_b_pts = b_rad_pts * 180/np.pi
            d_grid, l_grid, b_grid = np.meshgrid(self.density_grid_d_pts, self.density_grid_l_pts, 
                                                 self.density_grid_b_pts)
            grid_shape = d_grid.shape
            r_flat, phi_flat, z_flat = self.coord_trans.dlb_to_rphiz(d_grid.ravel(), 
                                                l_grid.ravel(), b_grid.ravel())

            # Evaluate the density across the grid, then integrate
            self.density_grid = self.density(r_flat, phi_flat, z_flat).reshape(grid_shape)
            vol_elem = d_grid**2 * np.cos(b_grid*np.pi/180)
            self.vol_elem=vol_elem
            self.density_int_b = simpson(self.density_grid*vol_elem, x=self.density_grid_b_pts*np.pi/180, axis=2)
            self.density_int_l = simpson(self.density_int_b, x=self.density_grid_l_pts*np.pi/180, axis=0)
            self.total_mass = simpson(self.density_int_l, x=self.density_grid_d_pts, axis=0)

    def draw_random_positions(self, n_stars: int = 1) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Draw points from the density distribution. Uses a grid of density points and
        CDF inversion in 3 dimensions.

        Parameters
        ----------
        dist_max : float [kpc]
            upper distance limit
        n_stars : int, None, optional
            number of stars drawn
            if None return one position as float

        Returns
        -------

        x : float, ndarray [kpc]
            Cartesian X coordinate (centered at the galactic center) of the drawn positions
        y : float, ndarray [kpc]
            Cartesian Y coordinate (centered at the galactic center) of the drawn positions
        z : float, ndarray [kpc]
            Cartesian Z coordinate (centered at the galactic center) of the drawn positions

        d_kpc : float, ndarray [kpc]
            distances of the drawn positions
        star_l_deg : float, ndarray [deg]
            galactic longitude of the drawn positions
        star_b_deg : float, ndarray [deg]
            galactic latitude of the drawn positions
        """
        if n_stars==0:
            return np.empty(0),np.empty(0),np.empty(0),np.empty(0),np.empty(0),np.empty(0)

        if self.field_shape == 'circle':
            # Select distances
            d_cum_dens = np.cumsum(self.density_int_st_dir)-self.density_int_st_dir[0]
            rand_pts = np.random.rand(n_stars)*d_cum_dens[-1]
            d_kpc = np.interp(rand_pts, d_cum_dens, self.density_grid_d_pts)

            # Select radial direction from center
            idx_d_nearest = np.argmin(np.abs(self.density_grid_d_pts - d_kpc[:, None]), axis=1)
            int_st_rad_idx = self.density_int_st_rad[:,idx_d_nearest]
            st_dir_cum_dens = (np.cumsum(int_st_rad_idx, axis=0) - int_st_rad_idx[0])
            if np.any(st_dir_cum_dens[-1]==0.0):
                # Deal with edge case
                bad_pts = np.where(st_dir_cum_dens[-1]==0.0)
                new_idx_d_nearest = idx_d_nearest[bad_pts]+1
                idx_d_nearest[bad_pts] = new_idx_d_nearest
                new_int_st_rad_idx = self.density_int_st_rad[:,new_idx_d_nearest]
                st_dir_cum_dens[:,bad_pts] = (np.cumsum(new_int_st_rad_idx, axis=0) - new_int_st_rad_idx[0])[:,None,:]
                assert ~np.any(st_dir_cum_dens[-1]==0.0)
            rand_pts = np.random.random(n_stars)*st_dir_cum_dens[-1]
            near_pts_hi = np.minimum(np.maximum((st_dir_cum_dens < rand_pts).sum(axis=0),1),len(self.density_grid_st_dir)-1)
            near_pts_lo = np.maximum(near_pts_hi-1,0)
            near_pts_lo_dens = st_dir_cum_dens[near_pts_lo,range(n_stars)]
            lin_fac = (rand_pts - near_pts_lo_dens) / (st_dir_cum_dens[near_pts_hi,range(n_stars)]-near_pts_lo_dens)
            assert ~np.any(np.isnan(lin_fac))
            if n_stars>0:
                assert np.all(rand_pts>=near_pts_lo_dens)
                assert np.all(rand_pts<=st_dir_cum_dens[near_pts_hi,range(n_stars)])

            star_st_dir = (1-lin_fac)*self.density_grid_st_dir[near_pts_lo] + lin_fac*self.density_grid_st_dir[near_pts_hi]

            # Select angular distance from center
            idx_st_dir_nearest = np.argmin(np.abs(self.density_grid_st_dir-star_st_dir[:, None]), axis=1)
            rho_grid_idx = self.density_grid[idx_st_dir_nearest,idx_d_nearest,:]
            st_rad_cum_dens = np.cumsum(rho_grid_idx, axis=1) - rho_grid_idx[:,0][:,None]
            if np.any(st_rad_cum_dens[-1]==0.0):
                # Deal with edge case
                bad_pts = np.where(st_rad_cum_dens[-1]==0.0)
                new_idx_st_dir_nearest = idx_st_dir_nearest[bad_pts]+1
                new_rho_grid_idx = self.density_grid[new_idx_st_dir_nearest,idx_d_nearest[bad_pts],:]
                st_rad_cum_dens[bad_pts] = (np.cumsum(new_rho_grid_idx, axis=1) - new_rho_grid_idx[:,0][:,None])[:,None,:]
                #print('did it work ?',~np.any(b_cum_dens[-1]==0.0))
            rand_pts = np.random.random(n_stars)*st_rad_cum_dens[:,-1]
            near_pts_hi = np.minimum(np.maximum((st_rad_cum_dens <= rand_pts[:,None]).sum(axis=1),1),len(self.density_grid_st_rad)-1)
            near_pts_lo = np.maximum(near_pts_hi-1,0)
            near_pts_lo_dens = st_rad_cum_dens[range(n_stars),near_pts_lo]
            lin_fac = np.nan_to_num((rand_pts - near_pts_lo_dens) / (st_rad_cum_dens[range(n_stars),near_pts_hi]-near_pts_lo_dens))
            assert np.all(lin_fac<=1)
            assert np.all(lin_fac>=0)
            star_st_rad = (1-lin_fac)*self.density_grid_st_rad[near_pts_lo] + lin_fac*self.density_grid_st_rad[near_pts_hi]

            # Get into physically meaningful coordinates
            delta_l_rad = star_st_rad * np.sin(star_st_dir)
            delta_b_rad = star_st_rad * np.cos(star_st_dir)
            l_rad, b_rad = self.rotate_00_to_lb(delta_l_rad, delta_b_rad)
            star_l_deg = l_rad * 180 / np.pi
            if np.abs(self.l_deg)<90:
                star_l_deg -= 360*(star_l_deg>180)
            star_b_deg = b_rad * 180 / np.pi

        elif self.field_shape == 'box':
            # Select distances
            d_cum_dens = np.cumsum(self.density_int_l)-self.density_int_l[0]
            rand_pts = np.random.rand(n_stars)*d_cum_dens[-1]
            d_kpc = np.interp(rand_pts, d_cum_dens, self.density_grid_d_pts)

            # Select longitudes
            idx_d_nearest = np.argmin(np.abs(self.density_grid_d_pts - d_kpc[:, None]), axis=1)
            int_b_idx = self.density_int_b[:,idx_d_nearest]
            l_cum_dens = (np.cumsum(int_b_idx, axis=0) - int_b_idx[0])
            if np.any(l_cum_dens[-1]==0.0):
                # Deal with edge case
                bad_pts = np.where(l_cum_dens[-1]==0.0)
                #print('fixing bad pts,', bad_pts, 'try upping dist idx by 1 for these')
                new_idx_d_nearest = idx_d_nearest[bad_pts]+1
                idx_d_nearest[bad_pts] = new_idx_d_nearest
                new_int_b_idx = self.density_int_b[:,new_idx_d_nearest]
                l_cum_dens[:,bad_pts] = (np.cumsum(new_int_b_idx, axis=0) - new_int_b_idx[0])[:,None,:]
                assert ~np.any(l_cum_dens[-1]==0.0)
            rand_pts = np.random.random(n_stars)*l_cum_dens[-1]
            near_pts_hi = np.minimum(np.maximum((l_cum_dens < rand_pts).sum(axis=0),1),len(self.density_grid_l_pts)-1)
            near_pts_lo = np.maximum(near_pts_hi-1,0)
            near_pts_lo_dens = l_cum_dens[near_pts_lo,range(n_stars)]
            lin_fac = (rand_pts - near_pts_lo_dens) / (l_cum_dens[near_pts_hi,range(n_stars)]-near_pts_lo_dens)
            assert ~np.any(np.isnan(lin_fac))
            if n_stars>0:
                assert np.all(rand_pts>=near_pts_lo_dens)
                assert np.all(rand_pts<=l_cum_dens[near_pts_hi,range(n_stars)])

            star_l_deg = (1-lin_fac)*self.density_grid_l_pts[near_pts_lo] + lin_fac*self.density_grid_l_pts[near_pts_hi]

            # Select latitudes
            idx_l_nearest = np.argmin(np.abs(self.density_grid_l_pts-star_l_deg[:, None]), axis=1)
            rho_grid_idx = self.density_grid[idx_l_nearest,idx_d_nearest,:]
            b_cum_dens = np.cumsum(rho_grid_idx, axis=1) - rho_grid_idx[:,0][:,None]
            if np.any(b_cum_dens[-1]==0.0):
                # Deal with edge case
                bad_pts = np.where(b_cum_dens[-1]==0.0)
                #pdb.set_trace()
                #print('fixing bad pts,', bad_pts, 'try upping l idx by 1 for these')
                new_idx_l_nearest = idx_l_nearest[bad_pts]+1
                new_rho_grid_idx = self.density_grid[new_idx_l_nearest,idx_d_nearest[bad_pts],:]
                b_cum_dens[bad_pts] = (np.cumsum(new_rho_grid_idx, axis=1) - new_rho_grid_idx[:,0][:,None])[:,None,:]
                #print('did it work ?',~np.any(b_cum_dens[-1]==0.0))
            rand_pts = np.random.random(n_stars)*b_cum_dens[:,-1]
            near_pts_hi = np.minimum(np.maximum((b_cum_dens <= rand_pts[:,None]).sum(axis=1),1),len(self.density_grid_b_pts)-1)
            near_pts_lo = np.maximum(near_pts_hi-1,0)
            near_pts_lo_dens = b_cum_dens[range(n_stars),near_pts_lo]
            lin_fac = np.nan_to_num((rand_pts - near_pts_lo_dens) / (b_cum_dens[range(n_stars),near_pts_hi]-near_pts_lo_dens))
            assert np.all(lin_fac<=1)
            assert np.all(lin_fac>=0)
            star_b_deg = (1-lin_fac)*self.density_grid_b_pts[near_pts_lo] + lin_fac*self.density_grid_b_pts[near_pts_hi]

        # estimate galactocentric coordinates
        x, y, z = self.coord_trans.dlb_to_xyz(d_kpc, star_l_deg, star_b_deg)

        return x, y, z, d_kpc, star_l_deg, star_b_deg

    def draw_random_point_in_slice(self, dist_inner: float, dist_outer: float, n_stars: int = 1,
                                    population_density_func=None) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        DEPRECATED / UNUSED IN VERSION 2
        Draw one or more random point in a slice
        given coordinates (self.l, self.b) [degrees],
        solid angle[steradians], and distance range[kpc]

        To get distance d, we draw from a cumulative quadratic distribution. To do so,
        we found the integrated r**2 such that our CDF is
        Prob(x) = (x**3 - dist_inner**3)/(dist_outer**3 - dist_inner**3)
        Then, we invert for x(Prob) so that we can draw Prob from Uniform(0,1)
        x = ((r_max**3 - r_min**3)*Prob + r_min**3)**(1/3)

        Parameters
        ----------
        dist_inner : float [kpc]
            lower distance
        dist_outer : float [kpc]
            upper distance
        n_stars : int, None, optional
            number of stars drawn
            if None return one position as float

        Returns
        -------

        x : float, ndarray [kpc]
            Cartesian X coordinate (centered at the galactic center) of the drawn positions
        y : float, ndarray [kpc]
            Cartesian Y coordinate (centered at the galactic center) of the drawn positions
        z : float, ndarray [kpc]
            Cartesian Z coordinate (centered at the galactic center) of the drawn positions

        d_kpc : float, ndarray [kpc]
            distances of the drawn positions
        star_l_deg : float, ndarray [deg]
            galactic longitude of the drawn positions
        star_b_deg : float, ndarray [deg]
            galactic latitude of the drawn positions
        """
        if n_stars==0:
            return np.empty(0),np.empty(0),np.empty(0),np.empty(0),np.empty(0),np.empty(0)

        # generate a cone uniformly around l=0, b=0 with solidAngle as covered area
        # Phi in the paper
        if (population_density_func is None) or (self.field_shape=='circle'):
            d_kpc = np.cbrt(np.random.uniform(dist_inner ** 3, dist_outer ** 3, size=n_stars))
            if self.field_shape=='circle':
                st_dir = np.random.uniform(0, 2 * np.pi, size=n_stars)
                # Theta in the paper
                st_rad = np.arccos(np.random.uniform(np.cos(self.lb_radius_deg*np.pi/180), 1, size=n_stars))
                # Estimate offset to center in ra and dec
                delta_l_rad = st_rad * np.sin(st_dir)
                delta_b_rad = st_rad * np.cos(st_dir)
            if self.field_shape=='box':
                delta_l_rad = np.pi/180 * self.l_length_deg/2 * np.random.uniform(-1, 1, size=n_stars)
                delta_b_rad = np.pi/180 * self.b_length_deg/2 * np.random.uniform(-1, 1, size=n_stars)

            star_l_rad, star_b_rad = self.rotate_00_to_lb(delta_l_rad, delta_b_rad)
            star_l_deg = star_l_rad * 180 / np.pi
            star_b_deg = star_b_rad * 180 / np.pi

        # estimate galactocentric coordinates
        x, y, z = self.coord_trans.dlb_to_xyz(d_kpc, star_l_deg, star_b_deg)

        return x, y, z, d_kpc, star_l_deg, star_b_deg

    def rotate_00_to_lb(self, delta_l: np.ndarray, delta_b: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray]:
        """
        Rotates coordinate system such that 0, 0 lands on self.l ,self.b

        Parameters
        ----------
        delta_l : float, ndarray [radians]
            difference in galactic longitude
        delta_b : float, ndarray [radians]
            difference in galactic longitude

        Returns
        -------
        star_l_rad : float, ndarray [radians]
            galactic longitude
        star_b_rad : float, ndarray [radians]
            galactic latitude
        """
        # calculate sin and cos.
        sin_theta, cos_theta = np.sin(delta_b), np.cos(delta_b)
        sin_phi, cos_phi = np.sin(delta_l), np.cos(delta_l)

        # estimate rotation matrix,
        # by rotating first around y-axis and then around z-axis
        mat = np.matmul(
            sp_utils.rotation_matrix(self.l_rad, axis='z'),
            sp_utils.rotation_matrix(self.b_rad, axis='y')
            )
        # convert to spherical coordinates
        vec = np.array([cos_theta * cos_phi, cos_theta * sin_phi, sin_theta])

        # apply rotation matrix
        vec = np.dot(mat, vec)

        # convert to galactic coordinates
        star_b_rad = np.arcsin(vec[2])
        star_l_rad = np.arctan2(vec[1], vec[0])
        star_l_rad += 2 * np.pi * (star_l_rad < 0)  # only if phi_prime_rad < 0
        # this way it works with both ndarray and floats

        return star_l_rad, star_b_rad
