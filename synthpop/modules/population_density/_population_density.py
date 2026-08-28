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
from scipy.integrate import cumulative_trapezoid, trapezoid

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
            coord_trans: ModuleType = sp_utils.coordinates_transformation.CoordTrans(),
            gamma_flare: float = None,
            h_flare: float = None,
            radius_flare: float = 0,
            grid_resolution: float = 0.010,
            logger: ModuleType = None,
            max_gc_dist: float = None,
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
        gamma_flare, h_flare, radius_flare: float
            parameters to implement the flare of the milky way
        grid_resolution: float
            required spatial resolution in kpc for the density integration/interpolation
            grid. for angular dimensions, use this cutoff at 10 kpc or max_distance if shorter.
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
        self.required_grid_resolution = grid_resolution
        self.initial_grid_resolution = np.maximum(0.01, grid_resolution)

        # Set limits for where we need to look for stars
        self.max_gc_dist = max_gc_dist

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
        Estimates the correction factor for the flare.
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
            if isinstance(field_scale_deg, (list, tuple, set)):
                self.l_length_deg = field_scale_deg[0]
                self.b_length_deg = field_scale_deg[1]
            else:
                self.l_length_deg = field_scale_deg
                self.b_length_deg = field_scale_deg
        self.max_distance = max_distance
        self.current_grid_resolution = self.initial_grid_resolution

        if self.logger is not None:
            self.logger.debug("setting up density grid")

        #Distance initial bounds-- not dependent on field shape
        d_min, d_max = 0, self.max_distance
        if self.max_gc_dist is not None:
            # Apply a 5% buffer to be safe
            d_min = self.sun.gal_dist-self.max_gc_dist*1.05
            d_max = self.sun.gal_dist+self.max_gc_dist*1.05

        if self.field_shape == 'circle':
            self.make_density_grid_circle(d_min, d_max, 0, 2*np.pi, 0, self.lb_radius_deg*np.pi/180)
            if self.total_mass>0.0:
                # Adjust distance grid for zero density regions if relevant
                nz_indices = np.flatnonzero(self.density_int_st_dir)
                dmin_idx = np.maximum(nz_indices[0]-1, 0)
                dmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_d_pts)-1)
                nz_indices = np.flatnonzero(trapezoid(self.density_int_st_rad, x=self.density_grid_d_pts, axis=1))
                stdir_min_idx = np.maximum(nz_indices[0]-1, 0)
                stdir_max_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_st_dir)-1)
                nz_indices = np.flatnonzero(trapezoid(trapezoid(self.density_grid_vscaled, x=self.density_grid_st_dir, 
                                                axis=0), x=self.density_grid_d_pts, axis=0))
                strad_min_idx = np.maximum(nz_indices[0]-1, 0)
                strad_max_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_st_rad)-1)
                while dmin_idx>0 or dmax_idx<(len(self.density_grid_d_pts)-1) \
                      or stdir_min_idx>0 or stdir_max_idx<(len(self.density_grid_st_dir)-1)  \
                      or strad_min_idx>0 or strad_max_idx<(len(self.density_grid_st_rad)-1)  :
                    self.current_grid_resolution =  np.maximum(self.current_grid_resolution/3,
                                                self.current_grid_resolution*(dmax_idx-dmin_idx)/len(self.density_grid_d_pts))
                    self.make_density_grid_circle(self.density_grid_d_pts[dmin_idx], self.density_grid_d_pts[dmax_idx],
                                                  self.density_grid_st_dir[stdir_min_idx], self.density_grid_st_dir[stdir_max_idx],
                                                  self.density_grid_st_rad[strad_min_idx], self.density_grid_st_rad[strad_max_idx])
                    nz_indices = np.flatnonzero(self.density_int_st_dir)
                    dmin_idx = np.maximum(nz_indices[0]-1, 0)
                    dmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_d_pts))
                    nz_indices = np.flatnonzero(trapezoid(self.density_int_st_rad, x=self.density_grid_d_pts, axis=1))
                    stdir_min_idx = np.maximum(nz_indices[0]-1, 0)
                    stdir_max_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_st_dir)-1)
                    nz_indices = np.flatnonzero(trapezoid(trapezoid(self.density_grid_vscaled, x=self.density_grid_st_dir, 
                                                    axis=0), x=self.density_grid_d_pts, axis=0))
                    strad_min_idx = np.maximum(nz_indices[0]-1, 0)
                    strad_max_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_st_rad)-1)
                    #pdb.set_trace()
                # If needed, zoom to required resolution
                if self.current_grid_resolution>self.required_grid_resolution:
                    self.current_grid_resolution = self.required_grid_resolution
                    self.make_density_grid_circle(self.density_grid_d_pts[dmin_idx], self.density_grid_d_pts[dmax_idx],
                                                  self.density_grid_st_dir[stdir_min_idx], self.density_grid_st_dir[stdir_max_idx],
                                                  self.density_grid_st_rad[strad_min_idx], self.density_grid_st_rad[strad_max_idx])

        elif self.field_shape == 'box':
            dl_min, dl_max = -self.l_length_deg/2 * np.pi/180, self.l_length_deg/2 * np.pi/180
            db_min, db_max = -self.b_length_deg/2 * np.pi/180, self.b_length_deg/2 * np.pi/180
            # Create the initial density grid
            self.make_density_grid_box(d_min, d_max, dl_min, dl_max, db_min, db_max)
            if (self.total_mass>0.0):
                # Adaptively adjust distance grid for zero density regions if relevant
                nz_indices = np.flatnonzero(self.density_int_l)
                dmin_idx = np.maximum(nz_indices[0]-1, 0)
                dmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_d_pts)-1)
                nz_indices = np.flatnonzero(trapezoid(self.density_int_b, x=self.density_grid_d_pts, axis=1))
                lmin_idx = np.maximum(nz_indices[0]-1, 0)
                lmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_dl_pts)-1)
                nz_indices = np.flatnonzero(trapezoid(trapezoid(self.density_grid_vscaled, x=self.density_grid_dl_pts, 
                                                axis=0), x=self.density_grid_d_pts, axis=0))
                bmin_idx = np.maximum(nz_indices[0]-1, 0)
                bmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_db_pts)-1)
                while dmin_idx>0 or dmax_idx<(len(self.density_grid_d_pts)-1) \
                        or lmin_idx>0 or lmax_idx<(len(self.density_grid_dl_pts)-1) \
                        or bmin_idx>0 or bmax_idx<(len(self.density_grid_db_pts)-1):
                    self.current_grid_resolution =  np.maximum(self.current_grid_resolution/3,
                                    self.current_grid_resolution*(dmax_idx-dmin_idx)/len(self.density_grid_d_pts))
                    self.make_density_grid_box(self.density_grid_d_pts[dmin_idx], self.density_grid_d_pts[dmax_idx],
                                               self.density_grid_dl_pts[lmin_idx], self.density_grid_dl_pts[lmax_idx],
                                               self.density_grid_db_pts[bmin_idx], self.density_grid_db_pts[bmax_idx])
                    nz_indices = np.flatnonzero(self.density_int_l)
                    dmin_idx = np.maximum(nz_indices[0]-1, 0)
                    dmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_d_pts)-1)
                    nz_indices = np.flatnonzero(trapezoid(self.density_int_b, x=self.density_grid_d_pts, axis=1))
                    lmin_idx = np.maximum(nz_indices[0]-1, 0)
                    lmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_dl_pts)-1)
                    nz_indices = np.flatnonzero(trapezoid(trapezoid(self.density_grid_vscaled, x=self.density_grid_dl_pts, 
                                                    axis=0), x=self.density_grid_d_pts, axis=0))
                    bmin_idx = np.maximum(nz_indices[0]-1, 0)
                    bmax_idx = np.minimum(nz_indices[-1]+1,len(self.density_grid_db_pts)-1)
                # If needed, get into required resolution
                if self.current_grid_resolution>self.required_grid_resolution:
                    self.current_grid_resolution = self.required_grid_resolution
                    self.make_density_grid_box(self.density_grid_d_pts[dmin_idx], self.density_grid_d_pts[dmax_idx],
                                               self.density_grid_dl_pts[lmin_idx], self.density_grid_dl_pts[lmax_idx],
                                               self.density_grid_db_pts[bmin_idx], self.density_grid_db_pts[bmax_idx])

    def make_density_grid_circle(self, d_min, d_max, stdir_min, stdir_max, strad_min, strad_max):
        # Set resolution
        grid_d_n_pts = int(np.ceil((d_max-d_min)/self.current_grid_resolution)) + 1
        grid_d_n_pts = np.maximum(grid_d_n_pts, 51)
        self.density_grid_d_pts = np.linspace(d_min, d_max, grid_d_n_pts)

        # Set up a grid
        res_dist = np.minimum(self.sun.gal_dist,d_max)
        grid_st_dir_n_pts = int(np.ceil(res_dist*strad_max/self.current_grid_resolution))
        grid_st_dir_n_pts = np.maximum(grid_st_dir_n_pts, 50)
        self.density_grid_st_dir = np.linspace(stdir_min, stdir_max, grid_st_dir_n_pts)
        grid_st_rad_n_pts = int(np.ceil(res_dist*(strad_max-strad_min)/self.current_grid_resolution))
        grid_st_rad_n_pts = np.maximum(grid_st_rad_n_pts,25)
        self.density_grid_st_rad = np.linspace(np.sqrt(strad_min), np.sqrt(strad_max), grid_st_rad_n_pts)**2
        d_grid, st_dir_grid, st_rad_grid = np.meshgrid(self.density_grid_d_pts, self.density_grid_st_dir, 
                                             self.density_grid_st_rad)
        self.current_grid_resolution = np.max([(self.density_grid_d_pts[1]-self.density_grid_d_pts[0]),
                                            (self.density_grid_st_dir[1]-self.density_grid_st_dir[0])*strad_max*res_dist,
                                            (self.density_grid_st_rad[-1]-self.density_grid_st_rad[-2])*res_dist])
        grid_shape = d_grid.shape

        # Get into physically meaningful coordinates
        delta_l_rad = st_rad_grid * np.sin(st_dir_grid)
        delta_b_rad = st_rad_grid * np.cos(st_dir_grid)
        l_grid_ravel, b_grid_ravel = self.rotate_00_to_lb(delta_l_rad.ravel(), delta_b_rad.ravel())
        r_flat, phi_flat, z_flat = self.coord_trans.dlb_to_rphiz(d_grid.ravel(), 
                                            l_grid_ravel*180/np.pi, b_grid_ravel*180/np.pi)

        # Get density at points and integrate
        self.density_grid = self.density(r_flat, phi_flat, z_flat).reshape(grid_shape)
        vol_elem = d_grid**2 * np.sin(st_rad_grid)
        self.density_grid_vscaled = self.density_grid*vol_elem
        self.density_int_st_rad = trapezoid(self.density_grid_vscaled, x=self.density_grid_st_rad, axis=2)
        self.density_int_st_dir = trapezoid(self.density_int_st_rad, x=self.density_grid_st_dir, axis=0)
        self.total_mass = trapezoid(self.density_int_st_dir, x=self.density_grid_d_pts, axis=0)

    def make_density_grid_box(self, d_min, d_max, dl_min, dl_max, db_min, db_max):
        # Set resolution
        grid_d_n_pts = int(np.ceil((d_max-d_min)/self.current_grid_resolution))+1
        grid_d_n_pts = np.maximum(grid_d_n_pts, 50)
        self.density_grid_d_pts = np.linspace(d_min, d_max, grid_d_n_pts)

        # Set up a grid
        res_dist = np.minimum(self.sun.gal_dist,d_max)
        delta_l_rad_n_pts = int(np.ceil(res_dist*(dl_max-dl_min)/self.current_grid_resolution))
        delta_l_rad_n_pts = np.maximum(40, delta_l_rad_n_pts)
        delta_b_rad_n_pts = int(np.ceil(res_dist*(db_max-db_min)/self.current_grid_resolution))
        delta_b_rad_n_pts = np.maximum(40, delta_b_rad_n_pts)
        self.density_grid_dl_pts = np.linspace(dl_min, dl_max, delta_l_rad_n_pts)
        self.density_grid_db_pts = np.linspace(db_min, db_max, delta_b_rad_n_pts)
        d_grid, dl_grid, db_grid = np.meshgrid(self.density_grid_d_pts, self.density_grid_dl_pts, 
                                             self.density_grid_db_pts)
        self.current_grid_resolution = np.max([(self.density_grid_d_pts[1]-self.density_grid_d_pts[0]),
                                            (self.density_grid_dl_pts[1]-self.density_grid_dl_pts[0])*res_dist,
                                            (self.density_grid_db_pts[1]-self.density_grid_db_pts[0])*res_dist])
        grid_shape = d_grid.shape

        # Get into physically meaningful coordinates
        l_grid_ravel, b_grid_ravel = self.rotate_00_to_lb(dl_grid.ravel(), db_grid.ravel())
        r_flat, phi_flat, z_flat = self.coord_trans.dlb_to_rphiz(d_grid.ravel(), 
                                            l_grid_ravel*180/np.pi, b_grid_ravel*180/np.pi)

        # Evaluate the density across the grid, then integrate
        self.density_grid = self.density(r_flat, phi_flat, z_flat).reshape(grid_shape)
        vol_elem = d_grid**2 * np.cos(b_grid_ravel.reshape(grid_shape)*np.pi/180)
        self.density_grid_vscaled = self.density_grid*vol_elem
        self.density_int_b = trapezoid(self.density_grid_vscaled, x=self.density_grid_db_pts, axis=2)
        self.density_int_l = trapezoid(self.density_int_b, x=self.density_grid_dl_pts, axis=0)
        self.total_mass = trapezoid(self.density_int_l, x=self.density_grid_d_pts, axis=0)

    def cumulative_integral(self, y, x=None, axis=-1, initial=0.0):
        # I had initially used simpson here, but it was unstable for the NSC
        # I think a more reliable method is to increase grid points and do trapezoid
        res = cumulative_trapezoid(y, x=x, axis=axis, initial=initial)
        return res

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
        if n_stars == 0:
            return np.empty(0), np.empty(0), np.empty(0), np.empty(0), np.empty(0), np.empty(0)

        idx_stars = np.arange(n_stars)
        if self.field_shape == 'circle':

            # First, select distances
            # Use quadratic CDF inversion w/ density point + cumulative density information
            #
            # Set up cumulative density grid and select surrounding points for each
            d_cum_dens = self.cumulative_integral(self.density_int_st_dir, x=self.density_grid_d_pts)
            rand_pts_d = np.random.rand(n_stars) * d_cum_dens[-1]
            near_d_hi = np.clip(np.searchsorted(d_cum_dens, rand_pts_d), 1, len(self.density_grid_d_pts) - 1)
            near_d_lo = near_d_hi - 1
            d1 = self.density_grid_d_pts[near_d_lo]
            d2 = self.density_grid_d_pts[near_d_hi]
            cum_dens_d_lo = d_cum_dens[near_d_lo]
            cum_dens_d_hi = d_cum_dens[near_d_hi]
            rho1_d = self.density_int_st_dir[near_d_lo]
            rho2_d = self.density_int_st_dir[near_d_hi]
            cum_dens_d_diff = cum_dens_d_hi - cum_dens_d_lo
            assert np.all(cum_dens_d_diff > 0.0), f"Density CDF must monotonically increase."
            # Get the random point in cumulative density and invert
            u_cell_d = np.clip((rand_pts_d - cum_dens_d_lo) / cum_dens_d_diff, 0.0, 1.0)
            delta_rho_d = rho2_d - rho1_d
            flat_mask_d = delta_rho_d==0.0
            inside_sqrt_d = np.maximum(0.0, rho1_d**2 + u_cell_d * (rho2_d**2 - rho1_d**2))
            frac_quad_d = (np.sqrt(inside_sqrt_d) - rho1_d) / np.where(flat_mask_d, 1.0, delta_rho_d)
            frac_d = np.where(flat_mask_d, u_cell_d, frac_quad_d)
            # Get the distance placement within the cell
            d_kpc = d1 + frac_d * (d2 - d1)

            # Second, select radial directions
            # Use quadratic CDF inversion w/ density point + cumulative density information
            #
            # Interpolate the density slices for the surrounding d_kpc
            t_d = (d_kpc - d1) / (d2 - d1)
            int_rad_lo = self.density_int_st_rad[:, near_d_lo]
            int_rad_hi = self.density_int_st_rad[:, near_d_hi]
            int_rad_interp = (1.0 - t_d) * int_rad_lo + t_d * int_rad_hi  
            # Set up cumulative density grid and select surrounding points for each
            st_dir_cum_dens = self.cumulative_integral(int_rad_interp, x=self.density_grid_st_dir, axis=0)
            rand_pts_dir = np.random.rand(n_stars) * st_dir_cum_dens[-1, :]
            near_dir_hi = np.clip((st_dir_cum_dens <= rand_pts_dir[np.newaxis, :]).sum(axis=0), 1, len(self.density_grid_st_dir) - 1)
            near_dir_lo = np.maximum(near_dir_hi - 1, 0)
            dens_near_dir_lo = self.density_grid_st_dir[near_dir_lo]
            dens_near_dir_hi = self.density_grid_st_dir[near_dir_hi]
            cum_dens_dir_lo = st_dir_cum_dens[near_dir_lo, idx_stars]
            cum_dens_dir_hi = st_dir_cum_dens[near_dir_hi, idx_stars]
            rho1_dir = int_rad_interp[near_dir_lo, idx_stars]
            rho2_dir = int_rad_interp[near_dir_hi, idx_stars]
            cum_dens_dir_diff = cum_dens_dir_hi - cum_dens_dir_lo
            assert np.all(cum_dens_dir_diff > 0.0), f"Directional CDF must monotonically increase."
            # Get the random point in cumulative density and invert
            u_cell_dir = np.clip((rand_pts_dir - cum_dens_dir_lo) / cum_dens_dir_diff, 0.0, 1.0)
            delta_rho_dir = rho2_dir - rho1_dir
            flat_mask_dir = delta_rho_dir==0.0
            inside_sqrt_dir = np.maximum(0.0, rho1_dir**2 + u_cell_dir * (rho2_dir**2 - rho1_dir**2))
            frac_quad_dir = (np.sqrt(inside_sqrt_dir) - rho1_dir) / np.where(flat_mask_dir, 1.0, delta_rho_dir)
            frac_dir = np.where(flat_mask_dir, u_cell_dir, frac_quad_dir)
            # Get the directional placement within the cell
            star_st_dir = dens_near_dir_lo + frac_dir * (dens_near_dir_hi - dens_near_dir_lo)

            # Third, we select angular radial distance
            # Use quadratic CDF inversion w/ density point + cumulative density information
            #
            # Interpolate accounting for star_st_dir and d_kpc
            idx_st_dir_hi = np.clip((self.density_grid_st_dir < star_st_dir[:, np.newaxis]).sum(axis=1), 1, len(self.density_grid_st_dir) - 1)
            idx_st_dir_lo = np.maximum(idx_st_dir_hi - 1, 0)
            st_dir_lo = self.density_grid_st_dir[idx_st_dir_lo]
            st_dir_hi = self.density_grid_st_dir[idx_st_dir_hi]
            w_dir = np.nan_to_num((star_st_dir - st_dir_lo) / (st_dir_hi - st_dir_lo))[:, np.newaxis]
            w_d = t_d[:, np.newaxis]
            rho_00 = self.density_grid_vscaled[idx_st_dir_lo, near_d_lo, :]
            rho_10 = self.density_grid_vscaled[idx_st_dir_hi, near_d_lo, :]
            rho_01 = self.density_grid_vscaled[idx_st_dir_lo, near_d_hi, :]
            rho_11 = self.density_grid_vscaled[idx_st_dir_hi, near_d_hi, :]
            rho_interp = ((1.0 - w_dir) * (1.0 - w_d) * rho_00 +
                           w_dir * (1.0 - w_d) * rho_10 +
                           (1.0 - w_dir) * w_d * rho_01 +
                           w_dir * w_d * rho_11)
            # Set up cumulative density grid and select surrounding points for each
            st_rad_cum_dens = self.cumulative_integral(rho_interp, x=self.density_grid_st_rad, axis=1)
            rand_pts_rad = np.random.rand(n_stars) * st_rad_cum_dens[:, -1]
            near_pts_hi = np.clip((st_rad_cum_dens <= rand_pts_rad[:, np.newaxis]).sum(axis=1), 1, len(self.density_grid_st_rad) - 1)
            near_pts_lo = np.maximum(near_pts_hi - 1, 0)
            r1 = self.density_grid_st_rad[near_pts_lo]
            r2 = self.density_grid_st_rad[near_pts_hi]
            cum_dens_rad_lo = st_rad_cum_dens[idx_stars, near_pts_lo]
            cum_dens_rad_hi = st_rad_cum_dens[idx_stars, near_pts_hi]
            rho1_rad = rho_interp[idx_stars, near_pts_lo]
            rho2_rad = rho_interp[idx_stars, near_pts_hi]
            cum_dens_rad_diff = cum_dens_rad_hi - cum_dens_rad_lo
            assert np.all(cum_dens_rad_diff > 0.0), f"Radial CDF must monotonically increase."
            # Get the random point in cumulative density and invert
            u_cell_rad = np.clip((rand_pts_rad - cum_dens_rad_lo) / cum_dens_rad_diff, 0.0, 1.0)
            delta_rho_rad = rho2_rad - rho1_rad
            flat_mask_rad = delta_rho_rad==0.0
            inside_sqrt_rad = np.maximum(0.0, rho1_rad**2 + u_cell_rad * (rho2_rad**2 - rho1_rad**2))
            frac_quad_rad = (np.sqrt(inside_sqrt_rad) - rho1_rad) / np.where(flat_mask_rad, 1.0, delta_rho_rad)
            frac_rad = np.where(flat_mask_rad, u_cell_rad, frac_quad_rad)
            # Get the angular radial placement within the cell
            star_st_rad = r1 + frac_rad * (r2 - r1)

            # Put into physical coordinates
            delta_l_rad = star_st_rad * np.sin(star_st_dir)
            delta_b_rad = star_st_rad * np.cos(star_st_dir)

        elif self.field_shape == 'box':

            # First, select distances
            # Use quadratic CDF inversion w/ density point + cumulative density information
            #
            # Set up cumulative density grid and select surrounding points for each
            d_cum_dens = self.cumulative_integral(self.density_int_l, x=self.density_grid_d_pts)
            rand_pts_d = np.random.rand(n_stars) * d_cum_dens[-1]
            near_d_hi = np.clip(np.searchsorted(d_cum_dens, rand_pts_d), 1, len(self.density_grid_d_pts) - 1)
            near_d_lo = near_d_hi - 1
            d1 = self.density_grid_d_pts[near_d_lo]
            d2 = self.density_grid_d_pts[near_d_hi]
            cum_dens_d_lo = d_cum_dens[near_d_lo]
            cum_dens_d_hi = d_cum_dens[near_d_hi]
            rho1_d = self.density_int_l[near_d_lo]
            rho2_d = self.density_int_l[near_d_hi]
            cum_dens_d_diff = cum_dens_d_hi - cum_dens_d_lo
            assert np.all(cum_dens_d_diff > 0.0), f"Density CDF must monotonically increase."
            # Get the random point in cumulative density and invert
            u_cell_d = np.clip((rand_pts_d - cum_dens_d_lo) / cum_dens_d_diff, 0.0, 1.0)
            delta_rho_d = rho2_d - rho1_d
            flat_mask_d = delta_rho_d==0
            inside_sqrt_d = np.maximum(0.0, rho1_d**2 + u_cell_d * (rho2_d**2 - rho1_d**2))
            frac_quad_d = (np.sqrt(inside_sqrt_d) - rho1_d) / np.where(flat_mask_d, 1.0, delta_rho_d)
            frac_d = np.where(flat_mask_d, u_cell_d, frac_quad_d)
            # Get the distance placement within the cell
            d_kpc = d1 + frac_d * (d2 - d1)

            # Second, select longitude offset
            # Use quadratic CDF inversion w/ density point + cumulative density information
            #
            # Interpolate the density slices for the surrounding d_kpc
            t_d = (d_kpc - d1) / (d2 - d1)
            int_l_lo = self.density_int_b[:, near_d_lo]
            int_l_hi = self.density_int_b[:, near_d_hi]
            int_l_interp = (1.0 - t_d) * int_l_lo + t_d * int_l_hi  
            # Set up cumulative density grid and select surrounding points for each
            l_cum_dens = self.cumulative_integral(int_l_interp, x=self.density_grid_dl_pts, axis=0)
            rand_pts_l = np.random.rand(n_stars) * l_cum_dens[-1, :]
            near_l_hi = np.clip((l_cum_dens <= rand_pts_l[np.newaxis, :]).sum(axis=0), 1, len(self.density_grid_dl_pts) - 1)
            near_l_lo = np.maximum(near_l_hi - 1, 0)
            l1 = self.density_grid_dl_pts[near_l_lo]
            l2 = self.density_grid_dl_pts[near_l_hi]
            cum_dens_l_lo = l_cum_dens[near_l_lo, idx_stars]
            cum_dens_l_hi = l_cum_dens[near_l_hi, idx_stars]
            rho1_l = int_l_interp[near_l_lo, idx_stars]
            rho2_l = int_l_interp[near_l_hi, idx_stars]
            cum_dens_l_diff = cum_dens_l_hi - cum_dens_l_lo
            assert np.all(cum_dens_l_diff > 0.0), f"Longitude CDF must monotonically increase."
            # Get the random point in cumulative density and invert
            u_cell_l = np.clip((rand_pts_l - cum_dens_l_lo) / cum_dens_l_diff, 0.0, 1.0)
            delta_rho_l = rho2_l - rho1_l
            flat_mask_l = delta_rho_l==0
            inside_sqrt_l = np.maximum(0.0, rho1_l**2 + u_cell_l * (rho2_l**2 - rho1_l**2))
            frac_quad_l = (np.sqrt(inside_sqrt_l) - rho1_l) / np.where(flat_mask_l, 1.0, delta_rho_l)
            frac_l = np.where(flat_mask_l, u_cell_l, frac_quad_l)
            # Get the longitude placement within the cell
            delta_l_rad = l1 + frac_l * (l2 - l1)

            # Third, select latitude offset
            # Use quadratic CDF inversion w/ density point + cumulative density information
            #
            # Interpolate accounting for delta_l_rad and d_kpc
            idx_l_hi = np.clip((self.density_grid_dl_pts < delta_l_rad[:, np.newaxis]).sum(axis=1), 1, len(self.density_grid_dl_pts) - 1)
            idx_l_lo = np.maximum(idx_l_hi - 1, 0)
            l_lo = self.density_grid_dl_pts[idx_l_lo]
            l_hi = self.density_grid_dl_pts[idx_l_hi]
            w_l = np.nan_to_num((delta_l_rad - l_lo) / (l_hi - l_lo))[:, np.newaxis]
            w_d = t_d[:, np.newaxis]
            rho_00 = self.density_grid_vscaled[idx_l_lo, near_d_lo, :]
            rho_10 = self.density_grid_vscaled[idx_l_hi, near_d_lo, :]
            rho_01 = self.density_grid_vscaled[idx_l_lo, near_d_hi, :]
            rho_11 = self.density_grid_vscaled[idx_l_hi, near_d_hi, :]
            rho_interp = ((1.0 - w_l) * (1.0 - w_d) * rho_00 +
                           w_l * (1.0 - w_d) * rho_10 +
                           (1.0 - w_l) * w_d * rho_01 +
                           w_l * w_d * rho_11)
            # Set up cumulative density grid and select surrounding points for each
            b_cum_dens = self.cumulative_integral(rho_interp, x=self.density_grid_db_pts, axis=1)

            rand_pts_b = np.random.rand(n_stars) * b_cum_dens[:, -1]
            near_b_hi = np.clip((b_cum_dens <= rand_pts_b[:, np.newaxis]).sum(axis=1), 1, len(self.density_grid_db_pts) - 1)
            near_b_lo = np.maximum(near_b_hi - 1, 0)
            b1 = self.density_grid_db_pts[near_b_lo]
            b2 = self.density_grid_db_pts[near_b_hi]
            cum_dens_b_lo = b_cum_dens[idx_stars, near_b_lo]
            cum_dens_b_hi = b_cum_dens[idx_stars, near_b_hi]
            rho1_b = rho_interp[idx_stars, near_b_lo]
            rho2_b = rho_interp[idx_stars, near_b_hi]
            cum_dens_b_diff = cum_dens_b_hi - cum_dens_b_lo
            assert np.all(cum_dens_b_diff > 0.0), f"Latitude CDF must monotonically increase."
            # Get the random point in cumulative density and invert
            u_cell_b = np.clip((rand_pts_b - cum_dens_b_lo) / cum_dens_b_diff, 0.0, 1.0)
            delta_rho_b = rho2_b - rho1_b
            flat_mask_b = delta_rho_b==0
            inside_sqrt_b = np.maximum(0.0, rho1_b**2 + u_cell_b * (rho2_b**2 - rho1_b**2))
            frac_quad_b = (np.sqrt(inside_sqrt_b) - rho1_b) / np.where(flat_mask_b, 1.0, delta_rho_b)
            frac_b = np.where(flat_mask_b, u_cell_b, frac_quad_b)
            # Get the latitude placement within the cell
            delta_b_rad = b1 + frac_b * (b2 - b1)

        # Convert from deltas from window center to real l & b
        star_l_rad, star_b_rad = self.rotate_00_to_lb(delta_l_rad, delta_b_rad)
        star_l_deg = star_l_rad * 180/np.pi
        star_b_deg = star_b_rad * 180/np.pi
        if np.abs(self.l_deg)<90:
            star_l_deg -= 360*(star_l_deg>180)
        # estimate galactocentric coordinates
        x, y, z = self.coord_trans.dlb_to_xyz(d_kpc, star_l_deg, star_b_deg)

        return x, y, z, d_kpc, star_l_deg, star_b_deg


    def draw_random_positions_rejection_method(self, n_stars: int = 1) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        FOR TESTING PURPOSES ONLY
        Draw points randomly uniformly then use rejection method to select according
        to density distribution.

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
        density_grid_max = np.max(self.density_grid)

        # Draw a bunch of points uniformly in the full window
        coords = np.array(self.draw_random_point_in_slice(0,self.max_distance, n_stars*1000))
        # Calculate density at these locations
        r, phi, z = self.coord_trans.dlb_to_rphiz(*coords[3:])
        rho_at_draws = self.density(r, phi, z)/density_grid_max
        # Draw random numbers to compare to the drawn densities to select kept stars
        keep_stars = np.where(rho_at_draws > np.random.rand(len(rho_at_draws)))[0]
        if len(keep_stars)>n_stars:
            keep_stars = keep_stars[:n_stars]
        # See if we have any stars to keep
        if len(keep_stars)==0:
            pos_list = np.array([[],[],[],[],[],[]])
        else:
            pos_list = coords[:, keep_stars]
        while len(pos_list[0]) < n_stars:
            #print(f"random positions drawn: {len(pos_list[0])} / {n_stars}")
            # Draw a bunch of points uniformly in the full window
            coords = np.array(self.draw_random_point_in_slice(0,self.max_distance, n_stars*100))
            # Calculate density at these locations
            r, phi, z = self.coord_trans.dlb_to_rphiz(*coords[3:])
            rho_at_draws = self.density(r, phi, z)/density_grid_max
            # Draw random numbers to compare to the drawn densities to select kept stars
            keep_stars = np.where(rho_at_draws > np.random.rand(len(rho_at_draws)))[0]
            if len(keep_stars)>(n_stars-len(pos_list[0])):
                keep_stars = keep_stars[:n_stars]
            if len(keep_stars)>0:
                pos_list = np.concatenate([pos_list, coords[:, keep_stars]], axis=1)
        return pos_list



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
