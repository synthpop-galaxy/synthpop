"""
Assign final compact object types and masses based on PopSyCLE (Rose et al 2022).
"""

__all__ = ["Spera15", ]
__author__ = "M.J. Huston"
__date__ = "2025-10-14"

import pandas
import numpy as np
from ._initial_final_mass_relation import InitialFinalMassRelation
#from scipy.stats import maxwell, uniform_direction
#from synthpop.synthpop_utils.coordinates_transformation import CoordTrans
import pdb
from typing import Set, Tuple, Dict, Union

class Spera15(InitialFinalMassRelation):
    """
    Post-processing to account for dim compact objects, based on the Spera15 IFMR from PopSyCLE (Rose et al 2022).
    """

    def __init__(self, logger, ifmr_name='SukhboldN20', **kwargs):
        super().__init__(logger, **kwargs)
        #: initial-final mass relation name to determine compact object masses.
        self.name='Spera15'
        self.spisea_ifmr_name = "IFMR_Spera15"

    def mass_spera15(self, m_zams, feh):
        """
        Remnant mass calculation from Spera et al 2015, appendix C
        Takes in the m_zams and Fe/H as lists
        """
        # First, calculate M_CO, based on M_ZAMS and Z
        # Note: z equation from Rose et al 2022,
        #       with z_sun from Ekstroem et al 2012
        z = 0.014*10**feh

        # First, calculate m_co

        # C13, 14, 15
        z_cat_1 = np.array([(z>4.0e-3).astype(int), 
                            ((z<=4.0e-3)&(z>=1.0e-3)).astype(int), 
                            ((z<1.0e-3)).astype(int)])
        b1s = np.array([59.63    - 2.969e3*z  + 4.988e4*z**2,     
                        40.98    + 3.415e4*z - 8.064e6*z**2,    
                        np.repeat(67.07,    len(z))])
        k1s = np.array([45.04    - 2.176e3*z  + 3.806e4*z**2,     
                        35.17    + 1.548e4*z - 3.759e6*z**2,    
                        np.repeat(46.89,    len(z))])
        k2s = np.array([138.9    - 4.664e3*z  + 5.106e4*z**2,     
                        20.36    + 1.162e5*z - 2.276e7*z**2,    
                        np.repeat(113.8,    len(z))])
        d1s = np.array([2.790e-2 - 1.780e-2*z + 77.05*z**2,       
                        2.500e-2 - 4.346*z   + 1.340e3*z**2,    
                        np.repeat(2.199e-2, len(z))])
        d2s = np.array([6.730e-3 + 2.690*z    - 52.39*z**2,       
                        1.750e-2 + 11.39*z   - 2.902e3*z**2,    
                        np.repeat(2.602e-2, len(z))])

        b1 = b1s[0]*z_cat_1[0] + b1s[1]*z_cat_1[1] + b1s[2]*z_cat_1[2]
        k1 = k1s[0]*z_cat_1[0] + k1s[1]*z_cat_1[1] + k1s[2]*z_cat_1[2]
        k2 = k2s[0]*z_cat_1[0] + k2s[1]*z_cat_1[1] + k2s[2]*z_cat_1[2]
        d1 = d1s[0]*z_cat_1[0] + d1s[1]*z_cat_1[1] + d1s[2]*z_cat_1[2]
        d2 = d2s[0]*z_cat_1[0] + d2s[1]*z_cat_1[1] + d2s[2]*z_cat_1[2]

        # C12
        g1 = 0.5 / (1 + 10**((k1-m_zams)*d1))
        g2 = 0.5 / (1 + 10**((k2-m_zams)*d2))
        # C11
        m_co = -2.0 + (b1+2.0)*(g1+g2)

        # Then, m_rem
        # Outer z condition for C1-3 vs C4-10
        z_cat_2 = np.array([(z<=5e-4).astype(int), (z>5e-4).astype(int)])
        # M_co condition for C4 and C1
        z_cat_2_1 = np.array([(m_co<5).astype(int), ((m_co>=5) & (m_co<10)).astype(int), (m_co>=10).astype(int)])
        # Inner z condition for C8-10
        z_cat_2_2 = np.array([(z>2e-3).astype(int), ((z<=2e-3) & (z>1e-3)).astype(int), (z<=1e-2).astype(int)])
        # Inner z condition for C6-7
        z_cat_2_3 = np.array([(z>1e-3).astype(int), (z<=1e-3).astype(int)])

        # m_rem for z<5e-4
        # C2-3
        m = -6.476e2*z + 1.911
        q = 2.300e3*z + 11.67
        p = -2.333 + 0.1559*m_co + 0.2700*m_co**2
        f = m*m_co + q
        m_rem_low_z = z_cat_2_1[0] * np.maximum(p, 1.27) + \
                      z_cat_2_1[1] * p + \
                      z_cat_2_1[2] * np.minimum(p, f)

        # m_rem for z>=5e-4
        m = z_cat_2_2[0]*np.repeat(1.217, len(z)) + z_cat_2_2[1]*(-43.82*z + 1.340)   + z_cat_2_2[2]*(-6.476e2*z + 1.911)
        q = z_cat_2_2[0]*np.repeat(1.061, len(z)) + z_cat_2_2[1]*(-1.296e4*z + 26.98) + z_cat_2_2[2]*(2.300e3*z + 11.67)
        a1  = z_cat_2_3[0]*(1.340 - 29.46 / (1 + (z/1.110e-3)**2.361))        + z_cat_2_3[1]*(1.105e5*z - 1.258e2)
        a2  = z_cat_2_3[0]*(80.22 - 74.73 * z**0.965 / (2.720e-3 + z**0.965)) + z_cat_2_3[1]*(91.56 - 1.957e4*z - 1.558e7*z**2)
        l   = z_cat_2_3[0]*(5.683 + 3.533 / (1 + (z/7.430e-3)**1.993))        + z_cat_2_3[1]*(1.134e4*z - 2.143)
        eta = z_cat_2_3[0]*(1.066 - 1.121 / (1 + (z/2.558e-2)**0.609))        + z_cat_2_3[1]*(3.090e-2 - 22.30*z + 7.363e4*z**2)

        h = a1 + (a2-a1)/(1+10**((l-m_co)*eta))
        f = m*m_co+q

        m_rem_high_z = z_cat_2_1[0] * np.maximum(h, 1.27) + \
                       z_cat_2_1[1] * h + \
                       z_cat_2_1[2] * np.maximum(h, f)

        m_rem = z_cat_2[0]*m_rem_low_z + z_cat_2[1]*m_rem_high_z

        return m_rem
            
    def compact_type_from_final(self, m_fin):
        """
        Determination of compact object type from final mass
        Based on PopSyCLE (Lam et al 2020; Rose et al 2022)
        Which pulls from Spera et al 2015

        Parameters
        ----------
        m_fin
            float value for final mass in units of solar mass
        Returns
        -------
        m_type
            integer value indicating object type
            1 = dim white dwarf
            2 = neutron star
            3 = black hole
        """
        return (m_fin<1.4).astype(int)*1 + ((m_fin>=1.4) & (m_fin<3)).astype(int) * 2 + (m_fin>=3).astype(int)*3

    def process_compact_objects(self, m_init: Union[np.ndarray, float],
                                      feh_init: Union[np.ndarray, float]):
        """
        Get the final masses and types for compact objects
        """
        # Cycle through evolved stars, calculating mass & type
        m_compact = self.mass_spera15(m_init,feh_init)
        m_type = self.compact_type_from_final(m_compact)
            
        return m_compact, 100+m_type
