""" Extinction maps from Surot et al 2020

"""

__all__ = ["Surot3d", ]
__author__ = "M.J. Huston"
__date__ = "2023-06-26"
__license__ = "GPLv3"
__version__ = "1.0.0"

import linecache
import numpy as np
import pandas as pd
from .. import const
from scipy.spatial import distance as sci_distance
from ._extinction import ExtinctionMap #, EXTINCTION_DIR
import time
import ebf
from scipy.interpolate import RegularGridInterpolator

class Surot3d(ExtinctionMap):
    """
    Extinction map from Surot et al. 2020
    Map must first be converted from E(J-Ks) to A_Ks values
    Note: 
        This is a quick implementation of a 2D extinction map which repurposes the methods desgined for a 3D map used in marshall.py. While the current version will primarily be used for extinction map comparison purposes, a more careful design will likely be beneficial in the future.  

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

    def __init__(self, dist_2d=8.15, **kwargs):
        # name of the extinction map used
        self.extinction_map_name = "Surot3d"
        # effective wavelength for VISTA K_s bandpass (http://svo2.cab.inta-csic.es/theory/fps/index.php?id=Paranal/VISTA.Ks&&mode=browse&gname=Paranal&gname2=VISTA#filter)
        self.ref_wavelength = 2.152152
        self.A_or_E_type = "A_Ks"
        self.dist_2d = dist_2d
        
        # placeholder for and location...
        self.l_deg = None
        self.b_deg = None
        self.extinction_in_map = None

        # get coordinates from surot_A_Ks_table1.csv
        tmp = pd.read_csv(f'{const.EXTINCTIONS_DIR}/surot_A_Ks_table_2D.csv',
            usecols=[0, 1, 2], sep=',', names=['l','b','A_Ks'])
        self.all_coords = np.transpose(np.array([tmp['l'],tmp['b']]))
        self.A_Ks_list = np.array(tmp['A_Ks'])

        self.l_stepsize = self.all_coords[1, 0] - self.all_coords[0, 0]
        
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

    def update_extinction_in_map(self, radius, force=False, **kwargs):
        """
        Check if the radius provided goes into a new bin, and if so update all the extinction values
        """
        # 3D portion for base fraction
        # Edge case: 0 for too low radius
        if self.l_deg<0:
            use_l = self.l_deg+360
        else:
            use_l = self.l_deg
        if radius < self.r_grid[0]:
            self.extinction_in_map = 0.0
        # Edge case: use extinction at furthest point for too high radius
        elif radius > self.r_grid[-1]:
            self.extinction_in_map = self.grid_interpolator_3d([use_l,self.b_deg, rgrid[-1]])[0]
        # Regular case, 3-D interpolation
        else:
            self.extinction_in_map = self.grid_interpolator_3d([use_l,self.b_deg, radius])[0]
            
        ext_norm = self.grid_interpolator_3d([use_l,self.b_deg, self.dist_2d])[0]
        
        self.extinction_in_map *= self.sightline[2]/ext_norm

        # a generic function we will place in all self.update_extinction functions
        return 1

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
        if self.l_deg>180:
            self.l_deg = 360.0-self.l_deg
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
        # find line of sight with the closest distance
        min_dist_arg = sci_distance.cdist(
                np.array([[self.l_deg, self.b_deg]]),
                self.all_coords
                ).argmin()

        self.sightline = [*self.all_coords[min_dist_arg], self.A_Ks_list[min_dist_arg]]
