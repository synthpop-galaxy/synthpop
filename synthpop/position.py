"""
This file includes the Position class.
It handles the generation of star positions within a given cone.
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
            circle radius or box half-widths
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
            self.l_hw_deg = None
            self.b_hw_deg = None
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
                self.l_hw_deg = field_scale_deg[0]
                self.b_hw_deg = field_scale_deg[1]
            else:
                self.l_hw_deg = field_scale_deg
                self.b_hw_deg = field_scale_deg

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

        # generate a cone uniformly around l=0, b=0 with solidAngle as covered area
        # Phi in the paper
        if True: #population_density_func is None:
            d_kpc = np.cbrt(np.random.uniform(dist_inner ** 3, dist_outer ** 3, size=n_stars))
            if self.field_shape=='circle':
                st_dir = np.random.uniform(0, 2 * np.pi, size=n_stars)
                # Theta in the paper
                st_rad = np.arccos(np.random.uniform(np.cos(self.lb_radius_deg*np.pi/180), 1, size=n_stars))
                # Estimate offset to center in ra and dec
                delta_l_rad = st_rad * np.sin(st_dir)
                delta_b_rad = st_rad * np.cos(st_dir)
            if self.field_shape=='box':
                delta_l_rad = np.pi/180 * self.l_hw_deg * np.random.uniform(-1, 1, size=n_stars)
                delta_b_rad = np.pi/180 * self.b_hw_deg * np.random.uniform(-1, 1, size=n_stars)
        else:
            pass
            #if self.field_shape == 'circle':
                # # Make a grid, sum over angles, inverse transform distances
                # dist_pts_unique = np.linspace(dist_inner, dist_outer, 20)
                # st_dir_pts_unique = np.linspace(0,2*np.pi, 20)
                # st_rad_pts_unique = np.linspace(0,self.lb_radius_deg, 20)
                # grid = np.meshgrid(dist_pts_unique, st_dir_pts_unique, st_rad_pts_unique)
                # dist_pts, st_dir_pts, st_rad_pts = np.fromiter(map(np.ravel, grid))
                # delta_l_pts_rad = st_rad_pts * np.sin(st_dir_pts)
                # delta_b_pts_rad = st_rad_pts * np.cos(st_dir_pts)
                # l_pts_rad, b_pts_rad = self.rotate_00_to_lb(delta_l_pts_rad, delta_b_pts_rad)
                # r_pts, phi_pts, z_pts = self.coord_trans.dlb_to_rphiz(dist_pts, l_pts_rad*180/np.pi, b_pts_rad*180/np.pi)
                # density_pts = population_density_func(r_pts, phi_pts, z_pts)
                # density_pts = density_pts.reshape(grid[0].shape)
                # density_dist_sum = np.cumsum(np.sum(density_pts, axis=(1,2)))
                # d_kpc = np.interp(np.random.uniform(0, density_dist_sum[-1], size=n_stars), density_dist_sum, dist_pts_unique)
                # # Then move on to the next dimension and do it again
                # grid = np.meshgrid(d_kpc, st_dir_pts_unique, st_rad_pts_unique)
                # dist_pts, st_dir_pts, st_rad_pts = np.fromiter(map(np.ravel, grid))
                # delta_l_pts_rad = st_rad_pts * np.sin(st_dir_pts)
                # delta_b_pts_rad = st_rad_pts * np.cos(st_dir_pts)
                # l_pts_rad, b_pts_rad = self.rotate_00_to_lb(delta_l_pts_rad, delta_b_pts_rad)
                # r_pts, phi_pts, z_pts = self.coord_trans.dlb_to_rphiz(dist_pts, l_pts_rad*180/np.pi, b_pts_rad*180/np.pi)
                # density_pts = population_density_func(r_pts, phi_pts, z_pts)
                # density_pts = density_pts.reshape(grid[0].shape)
                # density_st_dir_sum = np.cumsum(np.sum(density_pts, axis=2))
                # rand_pts = np.random.uniform(size=n_stars)
                # st_dir_radian = np.fromiter(map(np.interp(rand_pts[i]*density_str_dir_sum[-1, i], density_st_dir_sum[:,i], st_dir_pts_unique), range(n_stars)))
                # # Ok, one more dimension let's goooooo
                # grid1 = np.meshgrid(d_kpc, st_rad_pts_unique)
                # grid2 = np.meshgrid(st_dir_radian, st_rad_pts_unique)
                # dist_pts = np.ravel(grid1[0])
                # st_dir_pts = np.ravel(grid2[0])
                # st_rad_pts = np.ravel(grid2[1])
                # delta_l_pts_rad = st_rad_pts * np.sin(st_dir_pts)
                # delta_b_pts_rad = st_rad_pts * np.cos(st_dir_pts)
                # l_pts_rad, b_pts_rad = self.rotate_00_to_lb(delta_l_pts_rad, delta_b_pts_rad)
                # r_pts, phi_pts, z_pts = self.coord_trans.dlb_to_rphiz(dist_pts, l_pts_rad*180/np.pi, b_pts_rad*180/np.pi)
                # density_pts = population_density_func(r_pts, phi_pts, z_pts)
                # density_pts = density_pts.reshape(grid0[0].reshape)
                # density_st_rad_sum = np.cumsum(np.sum(density_pts))
                """Ok, so that was a mess. But it was close....."""

            #raise NotImplementedError('Density scaling within slice not yet implemented')
        # rotate cone to l_deg, b_deg
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
