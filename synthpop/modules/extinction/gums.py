"""
Extinction map used in the GAIA Universe Model, DR3 version.
Combines Lallement et al (2019) & Marshall (2006)

"""

__all__ = ["GUMS", ]
__author__ = "M.J. Huston"
__date__ = "2024-11-07"
__license__ = "GPLv3"
__version__ = "1.0.0"

import gzip
import h5py
import shutil
import numpy as np
from ._extinction import ExtinctionMap, EXTINCTION_DIR
from lallement import Lallement
from maps_from_dustmaps import MapsFromDustmaps
import time
import os
from .. import const
import requests

current_map_name = None
current_map_data = None

class GUMS(ExtinctionMap):
    """
    Extinction map from Lallement et al. 2019

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

    Methods
    -------

    update_extinction_in_map():
        placeholder for function that updates the total extinction or color excess
        in self.extinction_map_name

    get_map_properties():
        returns the basic parameters of the extinction map
        used for Communication between ExtinctionLaw and ExtinctionMap

    """

    def __init__(self, dr=0.001, return_functions=True, **kwargs):
        super().__init__(**kwargs)
        # name of the extinction map used
        self.extinction_map_name = "GUMS"
        self.return_functions = return_functions
        self.lallement = Lallement(dr=dr)
        self.marshall = MapsFromDustmaps(dustmap_name="marshall")
        self.ref_wavelength = self.lallement.ref_wavelength
        self.A_or_E_type = self.lallement.A_or_E_type
        # placeholder for and location...
        self.l_deg = None
        self.b_deg = None
        self.extinction_in_map = None
        
    def ext_func(self,l_deg,b_deg,dist):
        '''
        Get extinction value for multiple stars given their positions.
        '''
        # Draw a sightline outward for each star
        dist_max = np.max(dist)
        dist_pts = np.arange(0,dist_max,self.dr)[np.newaxis,:]
        l_rad, b_rad = l_deg[:,np.newaxis]*np.pi/180, b_deg[:,np.newaxis]*np.pi/180
        xm_pts = dist_pts*np.cos(b_rad)*np.cos(l_rad)
        ym_pts = dist_pts*np.cos(b_rad)*np.sin(l_rad)
        zm_pts = dist_pts*np.sin(b_rad)
        # Find where each line exits the grid
        within_grid = (np.abs(xm)<self.lallement.x_extent) & \
                        (np.abs(ym)<self.lallement.y_extent) & \
                        (np.abs(zm)<self.lallement.z_extent)
        end_pts = np.where(np.diff(within_grid))[1]
        end_dists = dist_pts[0][end_pts]
        # For stars within grid, use just Lallement.
        # For stars beyond, scale with addl distance as Marshall.
        value = self.lallement.ext_func(l_deg,b_deg,dist) * \ # Lallement value
                (self.marshall.get_map(l_deg,b_deg,dist) /  \ # Marshall value
                 self.marshall.get_map(l_deg,b_deg,end_dists)) ** \ # Marshall value @ grid edge
                (end_dists<dist) # Only apply Marshall if past grid
        return 

    def update_extinction_in_map(self, radius, force=False, **kwargs):
        """
        Returns the extinction for the current sight line and radial distance, or returns function to do so.

        Parameters
        ----------
        radius: float [kpc]
            radial distance of the current slice
        """

        if self.return_functions:
            self.extinction_in_map = self.ext_func
        else:
            self.extinction_in_map = self.ext_func(self.l_deg, self.b_deg, radius)
