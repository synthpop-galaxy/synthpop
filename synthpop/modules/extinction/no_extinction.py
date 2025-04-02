"""
Filler class to use for no extinciton
"""

__all__ = ["NoExtinction", ]
__author__ = "M.J. Huston"
__date__ = "2025-02-31"

import numpy as np
import pandas as pd
from ._extinction import ExtinctionMap

class NoExtinction(ExtinctionMap):
    """
    ExtinctionMap class to apply zero extinction in all directions
        
    Attributes
    ----------

    Methods
    -------
    extinction_in_map():
        function that returns 0 extinctions for all positions
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # name of the extinction map used
        self.extinction_map_name = "NoExtinction"
        self.ref_wavelength = 1.0
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
