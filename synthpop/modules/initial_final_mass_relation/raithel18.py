"""
Assign final compact object types and masses based on PopSyCLE (Rose et al 2022).
"""

__all__ = ["Raithel18", ]
__author__ = "M.J. Huston"
__date__ = "2025-10-14"

import pandas
import numpy as np
from ._initial_final_mass_relation import InitialFinalMassRelation
#from scipy.stats import maxwell, uniform_direction
#from synthpop.synthpop_utils.coordinates_transformation import CoordTrans
import pdb
from typing import Set, Tuple, Dict, Union

class Raithel18(InitialFinalMassRelation):
    """
    Post-processing to account for dim compact objects, based on Raithel18 PopSyCLE (Rose et al 2022).
    """

    def __init__(self, logger, ifmr_name='SukhboldN20', **kwargs):
        super().__init__(logger, **kwargs)
        #: initial-final mass relation name to determine compact object masses.
        self.name='Raithel18'
        self.spisea_ifmr_name = "IFMR_Raithel18"

    def mass_bh(self, m_zams, feh, f_ej=0.9):
        """
        Black hole mass calculation for Raithel18 and SukhboldN20

        Parameters
        ----------
        m_zams
            float or array of float values for initial stellar mass in units of solar mass
        f_ej
            float value representing the ejection fraction,
            or how much of the star's envelope is ejected in the supernova
            default value 0.9 adopted from Lam et al. (2020)

        Returns
        -------
        m_bh
            float or array of float values for final black hole mass in units of solar mass
        """
        m_bh_core_i = -2.049 + 0.4140 * m_zams
        m_bh_all_i = 15.52 - 0.3294 * (m_zams - 25.97) - 0.02121 * (
                    m_zams - 25.97) ** 2 + 0.003120 * (m_zams - 25.97) ** 3
        # branch ii
        m_bh_core_ii = 5.697 + 7.8598 * 10 ** 8 * m_zams ** -4.858
        # branch determination: 0 for i and 1 for ii
        branch = (m_zams > 42.21).astype(int)
        m_bh = (f_ej * m_bh_core_i + (1 - f_ej) * m_bh_all_i) * (1 - branch) + m_bh_core_ii * branch
        return m_bh

    def mass_ns(self, m_zams):
        """
        Neutron star final mass calculation, adopting the 1.36 Msun average
        with a standard deviation of 0.09.
        Based on PopSyCLE (Rose et al, 2022)

        Parameters
        ----------
        m_zams
            float or array of float values for initial stellar mass in units of solar mass

        Returns
        -------
        m_ns
            float or array of float values for final neutron star mass in units of solar mass
        """
        return np.random.normal(1.36, 0.09, len(m_zams))

    def mass_wd(self, m_zams):
        """
        White dwarf final mass calculation.
        Based on PopSyCLE (Rose et al. 2022)

        Parameters
        ----------
        m_zams
            float or array of float values for initial stellar mass in units of solar mass

        Returns
        -------
        m_wd
            float or array of float values for final white dwarf mass in units of solar mass
        """
        return 0.109 * m_zams + 0.394

    def compact_type_validation(self, m_type, m_prelim, m_ns):
        """
        Reassign low-mass BHs to NSs
        Based on PopSyCLE (Lam et al 2020; Rose et al 2022)

        Parameters
        ----------
        m_type
            int value for assigned compact object type
        m_prelim
            float value for compact object mass in units of solar mass
        m_ns
            float neutron star mass generated for the objects
        Returns
        -------
        m_type
            integer value indicating updated object type
            2 = neutron star
            3 = black hole
        m_compact
            updated mass of NS or BH
        """
        bh_to_ns = (m_type==3) & (m_prelim<3.0)
        m_type[bh_to_ns] = 2
        m_final = m_prelim
        m_final[bh_to_ns] = m_ns[bh_to_ns]
        return m_type, m_final

    def compact_type_from_initial(self, m_zams, feh):
        """
        Probabilistic drawing of compact object types
        Based on PopSyCLE (Lam et al 2020; Rose et al 2022)
        Which pulls from Rathiel et al 2018 and Sukhbold et al 2020

        Parameters
        ----------
        m_zams
            array of float values for initial stellar mass in units of solar mass
        feh
            array of float values for initial metallicity [Fe/H]
        Returns
        -------
        m_type
            array of integer values indicating object type
            0 = non-compact object or luminous white dwarf
            1 = dim white dwarf
            2 = neutron star
            3 = black hole
        """
        # Draw random numbers for bins that can be either NS or BH
        n_rand = np.random.uniform(size=len(m_zams))
        # Start with pre-CO objects, then go through mass bins, and assign appropriate type
        result = np.zeros(len(m_zams))
        result += ((m_zams>0.5)  & (m_zams<=9))    * 1
        result += ((m_zams>9)    & (m_zams<=15))   * 2
        result += ((m_zams>15)   & (m_zams<=17.8)) * ((n_rand<0.679)*2 + (n_rand>=0.679)*3)
        result += ((m_zams>17.8) & (m_zams<=18.5)) * ((n_rand<0.833)*2 + (n_rand>=0.833)*3)
        result += ((m_zams>18.5) & (m_zams<=21.7)) * ((n_rand<0.500)*2 + (n_rand>=0.500)*3)
        result += ((m_zams>21.7) & (m_zams<=25.2)) * 3
        result += ((m_zams>25.2) & (m_zams<=27.5)) * ((n_rand<0.652)*2 + (n_rand>=0.652)*3)
        result += ((m_zams>27.5) & (m_zams<=60))   * 3
        result += ((m_zams>60)   & (m_zams<=120))  * ((n_rand<0.400)*2 + (n_rand>=0.400)*3)
        return result

    def process_compact_objects(self, m_init: Union[np.ndarray, float],
                                      feh_init: Union[np.ndarray, float]):
        """
        Get the final masses and types for compact objects
        """
        # Probabilistic determination of object types
        m_type = self.compact_type_from_initial(m_init, feh_init)
        # Get possible masses and select by type
        m_wd = self.mass_wd(m_init)
        m_ns = self.mass_ns(m_init)
        m_bh = self.mass_bh(m_init, feh_init)
        m_compact = (m_bh * (m_type == 3).astype(int) +
                     m_ns * (m_type == 2).astype(int) +
                     m_wd * (m_type == 1).astype(int))
        m_type, m_compact = self.compact_type_validation(m_type, m_compact, m_ns)
            
        return m_compact, 100+m_type
