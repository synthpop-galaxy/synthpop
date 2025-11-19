"""
Binary companion generator, based on Raghavan et al. 2010

NOTE: The eccentricity is very hacky, there's not mathematical
ditstribution provided in the paper.
"""

__all__ = ["Multiplicity"]
__author__ = "M. Newman, M.J. Huston"
__credits__ = ["M. Newman, M.J. Huston"]
__date__ = "2025-10-15"

import pandas as pd
import numpy as np
from ._multiplicity import Multiplicity
import pdb
import synthpop.constants as const

try:
    from constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)
except (ImportError, ValueError):
    from synthpop.constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR, DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)

class Raghavan(Multiplicity):
    def __init__(self, **kwargs):
        """
        Hi
        """
        super().__init__(**kwargs)
        self.name='Raghavan'

    #@staticmethod
    def check_is_binary(self, pri_masses):
        """
        Probabilistic determination of binary star status based on temperature bins and probabilities.

        Parameters
        ----------
        pri_masses
            masses of the primaries

        Returns
        -------
        is_binary
            array of boolean values indicating whether each star is a binary
        """
        
        # Function from Duchene 2013
        binary_frac = (pri_masses<=0.1)*0.2 + \
                   (pri_masses>0.1)*0.3836*pri_masses**0.27

        is_binary = np.random.rand(len(pri_masses))>binary_frac

        return (is_binary, binary_frac)

    #@staticmethod
    def draw_companion_m_ratios(self, n):
        """
        Mass ratio (M2/M1) calculation from toy probability function based on Figure 16 of Raghavan et. al 2010

        Returns
        -------
        random_mass_ratio
            float value for the mass ratio (M2/M1) of the binary system
        """
        # Normalization factor (maximum value of N on Figure 16 of Raghavan 2010)
        #normalization_factor = 13

        # Separate systems into 3 different sections based on Raghavan 2010, Figure 16
        probability_bins = [0.103743, 0.778075+0.103743, 1.]    # Add previous bins to get *cumulative* values
        random_numbers = np.random.rand(n)
        bin_nos = np.searchsorted(probability_bins, random_numbers)

        random_mass_ratio = np.zeros(n)
        random_mass_ratio[bin_nos==0] = np.sqrt(2*np.random.rand(len(np.where(bin_nos==0)[0]))/50)
        random_mass_ratio[bin_nos==1] = np.random.uniform(0.2, 0.95, len(np.where(bin_nos==1)[0]))
        random_mass_ratio[bin_nos==2] = np.random.uniform(0.95, 1.0, len(np.where(bin_nos==2)[0]))

        return random_mass_ratio

    def draw_periods(self, n):
        """
        Find binary orbital period from Gaussian distribution in Figure 13 of Raghavan et al. 2010

        Returns
        -------
        periods
            float value of the period in days
        """
        # Draw a log(Period) from the Raghavan Figure 13 Gaussian
        logP = np.repeat(-1, n)
        mu = 5.03    # Value from Raghavan Figure 13
        sigma = 2.28    # Value from Raghavan Figure 13
        logP = np.random.normal(loc = mu, scale = sigma, size=n)

        # Keep drawing if we get a period less than 1 day (logP < 0)
        while np.any(logP <= 0) or np.any(logP >= 9):
            redraw_idx = ((logP <= 0) | (logP >= 9))
            logP[redraw_idx] = np.random.normal(loc=mu, scale=sigma, size=len(np.where(redraw_idx)[0]))

        return 10**logP

    def draw_eccentricities(self, periods):
        """
        Draw random eccentricities based on periods
        """
        e_arr = np.zeros(len(periods))
        not_circularized = (periods > 12)
        n_draw = np.sum(not_circularized)
        n_rand = np.random.rand(n_draw)
        e_drawn = n_rand * 0.6/0.8 * (n_rand<=0.8) + (1-np.sqrt((1-n_rand)/0.2)*0.4)*(n_rand>0.8)
        e_arr[not_circularized] = e_drawn
        return e_arr

    def generate_companions(self, pri_masses):
        """
        Generates companion stars

        Parameters
        ----------
        pri_masses : ndarray
            primary star masses

        Returns
        -------
        companions_table : dataframe
            companion properties
        """
                
        n_pri_stars = len(pri_masses)

        # Identify primary stars with companions
        (binary_flags, binary_frac) = self.check_is_binary(pri_masses)
        n_sec_stars = sum(binary_flags)
        pri_id = np.arange(n_pri_stars)[binary_flags]
        pri_masses_with_sec = pri_masses[binary_flags]
    
        # Draw an initial mass
        mass_ratios = self.draw_companion_m_ratios(n_sec_stars)
        sec_masses = mass_ratios * pri_masses_with_sec
        
        # Draw periods for binary stars
        periods = self.draw_periods(n_sec_stars)
        eccentricities = self.draw_eccentricities(periods)
                    
        return pri_id, sec_masses, periods, eccentricities
