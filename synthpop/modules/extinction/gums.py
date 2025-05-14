"""
Extinction map used in the GAIA Universe Model, DR3 version, which
combines Lallement et al (2019) & Marshall (2006)

Extinction is returned as A0, total extinction at reference wavelength 5500 Angstroms

GUMS documentation: https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu2sim/sec_cu2UM/ssec_cu2extmodel.html
"""

__all__ = ["Gums", ]
__author__ = "M.J. Huston"
__date__ = "2024-11-07"

import gzip
import h5py
import shutil
import numpy as np
try:
    from ._extinction import ExtinctionMap
    from .lallement import Lallement
    from .. import const
except ImportError:
    from _extinction import ExtinctionMap
    from lallement import Lallement
    import constants as const
import time
import os
import requests
import dustmaps.marshall
import astropy.units as u
from astropy.coordinates import SkyCoord

# empty dictionary to store dustmaps query
_query_dict = {}

class Gums(Lallement,ExtinctionMap):
    """
    Extinction map from Gaia Universe Model Snapshot

    Attributes
    ----------
    dr=0.001 : float
        step size for extinction integration in kpc
    per_sightline=True : boolean
        calculate extinction integral once per sightline if True (faster)
        calculate extinciton integral individually per star if False (slower)

    Methods
    -------
    gums_ext_func(l_deg, b_deg, dist):
        get extinction value in combo Lallement/Marshall combo map for list of star locations
    extinction_in_map(l_deg, b_deg, dist):
        equivalent to gums_ext_func
    """

    def __init__(self, dr=0.001, **kwargs):
        # Start from Lallement law
        super().__init__(dr=dr, **kwargs)
        self.extinction_map_name = "GUMS"
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
        """
        Estimates the extinction for a list of star positions.

        Parameters
        ----------
        l_deg: ndarray [degrees]
            galactic longitude
        b_deg: ndarray [degrees]
            galactic latitude
        dist: ndarray [kpc]
            radial distance from the Sun
        
        Returns
        -------
        extinction_value: ndarray [mag]
            extinction at each star position
        """
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
        return np.maximum(np.random.normal(value, value*0.1), 0)

    def extinction_in_map(self,l_deg,b_deg,dist):
        """
        Estimates the extinction for a list of star positions.

        Parameters
        ----------
        l_deg: ndarray [degrees]
            galactic longitude
        b_deg: ndarray [degrees]
            galactic latitude
        dist: ndarray [kpc]
            radial distance from the Sun
        
        Returns
        -------
        extinction_value: ndarray [mag]
            extinction at each star position
        """
        return self.gums_ext_func(l_deg, b_deg, dist)
