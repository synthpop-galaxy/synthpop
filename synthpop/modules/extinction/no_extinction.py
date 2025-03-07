"""
Filler class to use to not apply extinction
"""

__all__ = ["NoExtinction", ]
__author__ = "M.J. Huston"
__date__ = "2025-02-31"
__license__ = "GPLv3"
__version__ = "1.0.0"

import numpy as np
import pandas as pd
from ._extinction import ExtinctionMap

class NoExtinction(ExtinctionMap):
    """
    ExtinctionMap class to apply zero extinction in all directions
        
    Attributes
    ----------
    extinction_map_name : str
        name of the Extinction Map
    ref_wavelength : float
        reference wavelength for the extinction
    A_or_E_type : str
        Output type from the extinction map.
        If it starts with "A", A_or_E is handled  as a total extinction.
        If it starts with "E": A_or_E is handled as a color excess.

    Methods
    -------
    extinction_in_map():
        function that returns extinction at input positions
    get_map_properties():
        returns the basic parameters of the extinction map
        used for Communication between ExtinctionLaw and ExtinctionMap
    """

    def __init__(self, project_3d=True, dist_2d=8.15, **kwargs):
        super().__init__(**kwargs)
        # name of the extinction map used
        self.extinction_map_name = "NoExtinction"
        self.ref_wavelength = 2.152152
        self.A_or_E_type = "A_None"

    def extinction_in_map(self, l_deg, b_deg, dist):
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
        extinction at each star position defined as self.A_or_E_type
    """
        return dist*0.0
