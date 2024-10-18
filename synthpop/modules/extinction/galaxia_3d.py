""" 
Extinction from Galaxia, based on Schlegel et al 1998 2-D map, with 3-D dust disk model
NOTE: uses extinction_in_map value form; consider function form
NOTE: still undergoing testing

"""

__all__ = ["Galaxia3D", ]
__author__ = "M.J. Huston"
__date__ = "2024-04-18"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from .. import const
from ._extinction import ExtinctionMap #, EXTINCTION_DIR
import time
import ebf

class Galaxia3D(ExtinctionMap):
    """
    Extinction map used in Galaxia


    Attributes
    ----------
    extinction_map_name : str
        name of the Extinction Map
    l_deg : float
        galactic longitude in degree set by "update_sight-line"
    b_deg : float
        galactic latitude in degree set by "update_sight-line"

    ref_wavelength : float
        reference wavelength for the extinction

    A_or_E : float or function
        total extinction or color excess, from the extinction map.
        if it is a function it will be called

    A_or_E_type : str
        Output type from the extinction map.
        If it starts with "A", A_or_E is handled  as a total extinction.
        If it starts with "E": A_or_E is handled as a color excess.
    R_V : float
        interstellar reddening parameter

    Methods
    -------
    update_line_of_sight(l_deg, b_deg) :
        specified the new galactic coordinates
        calls find_sightline() and init_sightline()
    update_extinction_in_map() :

    init_sightline(self) :

    find_sightline(self) :

    update_extinction_in_map():
        placeholder for function that updates the total extinction or color excess
        in self.extinction_map_name
    set_R_V(self,R_V):
        set interstellar reddening parameter

    get_map_properties():
        returns the basic parameters of the extinction map
        used for Communication between ExtinctionLaw and ExtinctionMap

    """

    def __init__(self, which_2d='Shlegel',**kwargs):
        # name of the extinction map used
        self.extinction_map_name = "Galaxia3D"
        self.ref_wavelength = 0.4361
        self.ref_wavelength2 = 0.5448
        self.A_or_E_type = 'E(B-V)' 

        # placeholder for the coordinates of the sight-line
        # set by update_line_of_sight
        self.l_deg = None
        self.b_deg = None
        self.extinction_in_map = None
        self.base_extinction = None

        # Set up 3D grid
        mapfile_3d = ebf.read(f'{const.EXTINCTIONS_DIR}/Galaxia_ExMap3d_1024.ebf')
        map_grid_3d = mapfile_3d['exmap3d.xmms']
        map_data_3d = mapfile_3d['exmap3d.data']
        # Set up 3D interpolation
        l_grid_3d = np.append(np.arange(*map_grid_3d[0]),map_grid_3d[0][1])
        b_grid_3d = np.append(np.arange(*map_grid_3d[1]),map_grid_3d[1][1])
        self.r_grid = 10**np.append(np.arange(*map_grid_3d[2]),map_grid_3d[2][1])
        self.grid_interpolator_3d = RegularGridInterpolator((l_grid_3d,b_grid_3d,self.r_grid), 
            map_data_3d)

        # Set up 3d grid
        if which_2d=='Shlegel':
            mapfile_2d = ebf.read(f'{const.EXTINCTIONS_DIR}/Galaxia_Schlegel_4096.ebf')
        elif which_2d=='Solar':
            mapfile_2d = mapfile_3d
        map_grid_2d = mapfile_2d['exmap2d.xmms']
        map_data_2d = mapfile_2d['exmap2d.data']
        # set up interpolator
        l_grid_2d = np.append(np.arange(*map_grid_2d[0]),map_grid_2d[0][1])
        b_grid_2d = np.append(np.arange(*map_grid_2d[1]),map_grid_2d[1][1])
        self.grid_interpolator_2d = RegularGridInterpolator((l_grid_2d,b_grid_2d), 
            map_data_2d)
            
    def update_line_of_sight(self, l_deg, b_deg):
        """
        Set a new sight-line

        Parameters
        ----------
        l_deg : float [degree]
            galactic longitude
        b_deg : float [degree]
            galactic latitude

        """
        self.l_deg = l_deg
        if self.l_deg<0:
            self.l_deg+=360.0
        self.b_deg = b_deg
        self.find_sightline()
        self.init_sightline()
            
    def init_sightline(self):
        """
        Routine to prep the sightline, get first values, etc.
        """
        self.sightline_l_deg = self.sightline[0]
        self.sightline_b_deg = self.sightline[1]
        self.A_Ks = 0
        self.update_extinction_in_map(0)

    def find_sightline(self):
        """
        A function for the Surot map that finds the closest
        sightline to the coordinates specified.
        
        Can be used to reset the class to a new sightline,
        could be more efficient.
        """
        self.sightline = [self.l_deg, self.b_deg, self.grid_interpolator_2d([self.l_deg,self.b_deg])[0]]

    def update_extinction_in_map(self, radius: float):
        """
        Estimates the extinction for the current sight-line and radial distance
        store the result into self.extinction_in_map.

        Parameters
        ----------
        radius: float [kpc]
            radial distance of the current slice

        """

        # 3D portion for base fraction
        # Edge case: 0 for too low radius
        if radius < self.r_grid[0]:
            self.extinction_in_map = 0.0
        # Edge case: use extinction at furthest point for too high radius
        elif radius > self.r_grid[-1]:
            self.extinction_in_map = self.grid_interpolator_3d([self.l_deg,self.b_deg, rgrid[-1]])[0]
        # Regular case, 3-D interpolation
        else:
            self.extinction_in_map = self.grid_interpolator_3d([self.l_deg,self.b_deg, radius])[0]

        # 2D scaling
        self.extinction_in_map *= self.sightline[2]




