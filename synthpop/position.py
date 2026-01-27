"""
This file includes the Position class.
It handles the generation of star positions within a given window.
"""

__all__ = ['Position']
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-07-06"

from typing import Tuple
import numpy as np
try:
    from . import synthpop_utils as sp_utils
except ImportError:
    import synthpop_utils as sp_utils
import pdb
from scipy.integrate import simpson

class Position:
    """
    Position class for a Population class.
    This contains methods to randomly generate positions within a field/slice,

    Attributes
    ----------
    l_deg : float [degree]
        Galactic longitude of the current field in degrees.
    l_rad : float [radian]
        Galactic longitude of the current field in radians.
    b_deg : float [degree]
        Galactic latitude of the current field in degrees.
    b_rad : float [radian]
        Galactic latitude of the current field in radians.

    Methods
    -------
    __init__(l_deg: float, b_deg: float, solid_angle_sr: float) : None
        Initialize the Class
    update(*args, **kwargs) : None
        Update the class with new coordinates, passes arguments to __init__()
    draw_random_point_in_slice(dist_inner: float, dist_outer: float, N: int = 1): tuple
        generate N points within the slice
    rotate_00_to_lb(delta_l: ndarray, delta_b: ndarray) : Tuple[ndarray, ndarray]
        rotates a cone from pointing toward 0,0 to (l,b)
    """

    # placeholder for transformation matrix between the rectangular equatorial coordinate system
    # and the rectangular Galactic coordinate

    def __init__(self, logger, coord_trans, field_shape, field_scale_deg):
        """
        Initialization

        Parameters
        ----------
        logger
        coord_trans
        field_shape : str
            circle or box field shape
        field_scale_deg : float or np.ndarray
            circle radius or box side lengths
        """

        self.coord_trans = coord_trans
        # define coordinates of center of the field in degrees and radians
        self.l_deg = None
        self.l_rad = None
        self.b_deg = None
        self.b_rad = None
        self.logger = logger
        #  convert solid angle to half cone angle using wiki formula:
        self.field_shape = field_shape
        if field_shape=='circle':
            self.lb_radius_deg = None
        elif field_shape=='box':
            self.l_length_deg = None
            self.b_length_deg = None
        else:
            raise ValueError(f"field_shape {field_shape} not valid. Please use 'circle' or 'box'.")

    def update_location(self, l_deg: float, b_deg: float, field_shape: str, field_scale_deg: float):
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

    def draw_random_point_in_slice(self, dist_inner: float, dist_outer: float, n_stars: int = 1,
                                    population_density_func=None) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
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
        else:
            #if self.field_shape == 'circle'
            # need to get this method ready tooooooo
            if self.field_shape == 'box':
                # Box should be easier, right ?? RIGHT ???????
                # Let's make a grid
                d_pts = np.linspace(dist_inner, dist_outer, 501)
                l_pts = np.linspace(self.l_deg-self.l_length_deg/2/np.cos(self.b_rad), 
                                    self.l_deg+self.l_length_deg/2/np.cos(self.b_rad), 201)
                b_pts = np.linspace(self.b_deg-self.b_length_deg/2, self.b_deg+self.b_length_deg/2, 201)
                d_grid, l_grid, b_grid = np.meshgrid(d_pts, l_pts, b_pts)
                grid_shape = d_grid.shape
                d_flat = d_grid.ravel(); l_flat = l_grid.ravel(); b_flat = b_grid.ravel()
                r_flat, phi_flat, z_flat = self.coord_trans.dlb_to_rphiz(d_flat, l_flat, b_flat)
                r_grid = r_flat.reshape(grid_shape); phi_grid = phi_flat.reshape(grid_shape)
                z_grid = z_flat.reshape(grid_shape)
                # x_flat, y_flat, z_flat = self.coord_trans.dlb_to_xyz(d_flat, l_flat, b_flat)
                # x_grid = x_flat.reshape(grid_shape); y_grid = y_flat.reshape(grid_shape)

                # Evaluate the density across the grid, then integrate
                rho_grid = population_density_func(r_flat, phi_flat, z_flat).reshape(grid_shape)
                vol_elem = d_grid**2 * np.cos(b_grid*np.pi/180)
                int_b = simpson(rho_grid*vol_elem, x=b_pts*np.pi/180, axis=2)
                int_l = simpson(int_b, x=l_pts*np.pi/180, axis=0)
                int_d = simpson(int_l, x=d_pts, axis=0)
                # Dimensions aren't working  how I expected...... need to check for bugs

                # Select distances
                d_cum_dens = np.cumsum(int_l)-int_l[0]
                rand_pts = np.random.rand(n_stars)*d_cum_dens[-1]
                d_kpc = np.interp(rand_pts, d_cum_dens, d_pts)

                # Select longitudes
                idx_d_nearest = np.argmin(np.abs(d_pts - d_kpc[:, None]), axis=1)
                int_b_idx = int_b[:,idx_d_nearest]
                l_cum_dens = (np.cumsum(int_b_idx, axis=0) - int_b_idx[0])
                if np.any(l_cum_dens[-1]==0.0):
                    # Deal with edge case
                    bad_pts = np.where(l_cum_dens[-1]==0.0)
                    #print('fixing bad pts,', bad_pts, 'try upping dist idx by 1 for these')
                    new_idx_d_nearest = idx_d_nearest[bad_pts]+1
                    idx_d_nearest[bad_pts] = new_idx_d_nearest
                    new_int_b_idx = int_b[:,new_idx_d_nearest]
                    l_cum_dens[:,bad_pts] = (np.cumsum(new_int_b_idx, axis=0) - new_int_b_idx[0])[:,None,:]
                    assert ~np.any(l_cum_dens[-1]==0.0)
                rand_pts = np.random.random(n_stars)*l_cum_dens[-1]
                near_pts_hi = np.minimum(np.maximum((l_cum_dens < rand_pts).sum(axis=0),1),len(l_pts)-1)
                near_pts_lo = np.maximum(near_pts_hi-1,0)
                near_pts_lo_dens = l_cum_dens[near_pts_lo,range(n_stars)]
                lin_fac = (rand_pts - near_pts_lo_dens) / (l_cum_dens[near_pts_hi,range(n_stars)]-near_pts_lo_dens)
                assert ~np.any(np.isnan(lin_fac))
                if n_stars>0:
                    assert np.all(rand_pts>=near_pts_lo_dens)
                    assert np.all(rand_pts<=l_cum_dens[near_pts_hi,range(n_stars)])

                star_l_deg = (1-lin_fac)*l_pts[near_pts_lo] + lin_fac*l_pts[near_pts_hi]

                # Select latitudes
                idx_l_nearest = np.argmin(np.abs(l_pts-star_l_deg[:, None]), axis=1)
                rho_grid_idx = rho_grid[idx_l_nearest,idx_d_nearest,:]
                b_cum_dens = np.cumsum(rho_grid_idx, axis=1) - rho_grid_idx[0]
                if np.any(b_cum_dens[-1]==0.0):
                    # Deal with edge case
                    bad_pts = np.where(b_cum_dens[-1]==0.0)
                    #pdb.set_trace()
                    #print('fixing bad pts,', bad_pts, 'try upping l idx by 1 for these')
                    new_idx_l_nearest = idx_l_nearest[bad_pts]+1
                    new_rho_grid_idx = rho_grid[new_idx_l_nearest,idx_d_nearest[bad_pts],:]
                    b_cum_dens[bad_pts] = (np.cumsum(new_rho_grid_idx, axis=1) - new_rho_grid_idx[0])[:,None,:]
                    #print('did it work ?',~np.any(b_cum_dens[-1]==0.0))
                rand_pts = np.random.random(n_stars)*b_cum_dens[:,-1]
                near_pts_hi = np.minimum(np.maximum((b_cum_dens <= rand_pts[:,None]).sum(axis=1),1),len(b_pts)-1)
                near_pts_lo = np.maximum(near_pts_hi-1,0)
                near_pts_lo_dens = b_cum_dens[range(n_stars),near_pts_lo]
                lin_fac = np.nan_to_num((rand_pts - near_pts_lo_dens) / (b_cum_dens[range(n_stars),near_pts_hi]-near_pts_lo_dens))
                star_b_deg = (1-lin_fac)*b_pts[near_pts_lo] + lin_fac*b_pts[near_pts_hi]

                print('coords ready')

            #raise NotImplementedError('Density scaling within slice not yet implemented')

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
