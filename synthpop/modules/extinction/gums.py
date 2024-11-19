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
import time
import os
from .. import const
import requests
import dustmaps.marshall
import astropy.units as u
from astropy.coordinates import SkyCoord

# empty dictionary to store dustmaps query
_query_dict = {}

class Gums(Lallement,ExtinctionMap):
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
        # Start from Lallement law
        super().__init__(**kwargs)
        self.extinction_map_name = "GUMS"
        self.return_functions=return_functions
        # Set up Marshall map
        if not os.path.isdir(os.path.join(dustmaps.std_paths.data_dir(), "marshall")):
            url = 'https://dustmaps.readthedocs.io/en/latest/installation.html'
            module = dustmaps.marshall.MarshallQuery.__module__
            print("Downloading Marshall map")
            dustmaps.marshall.fetch()
        global _query_dict
        if 'marshall' not in _query_dict:
            _query_dict['marshall'] = dustmaps.marshall.MarshallQuery
        self.marshall_query = _query_dict['marshall']()
        
    def gums_ext_func(self,l_deg,b_deg,dist):
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
        within_grid = (np.abs(dist_pts*np.cos(b_rad)*np.cos(l_rad))<self.x_extent) & \
                        (np.abs(dist_pts*np.cos(b_rad)*np.sin(l_rad))<self.y_extent) & \
                        (np.abs(dist_pts*np.sin(b_rad))<self.z_extent)
        star_idx,end_idx= np.where(np.diff(within_grid))
        end_dists = 1.0*dist
        #print(star_idx,end_idx)
        #print(dist_pts)
        if len(star_idx>0):
            end_dists[star_idx] = dist_pts[0][end_idx]
        #print(end_dists)
        # For stars within grid, use just Lallement.
        # For stars beyond and in plane, scale with addl distance as Marshall.
        value = self.lallement_ext_func(l_deg,b_deg,dist) * \
                np.nan_to_num(self.marshall_query(SkyCoord(l_deg*u.deg,b_deg*u.deg,distance=dist*u.kpc, frame='galactic')) /  \
                 self.marshall_query(SkyCoord(l_deg*u.deg,b_deg*u.deg,distance=end_dists*u.kpc, frame='galactic')), nan=1.0) ** \
                (end_dists<dist)
        return np.maximum(value + np.random.normal(value, value*0.1), 0)

    def update_extinction_in_map(self, radius, force=False, **kwargs):
        """
        Returns the extinction for the current sight line and radial distance, or returns function to do so.

        Parameters
        ----------
        radius: float [kpc]
            radial distance of the current slice
        """

        if self.return_functions:
            self.extinction_in_map = self.gums_ext_func
        else:
            self.extinction_in_map = self.gums_ext_func(self.l_deg, self.b_deg, radius)
