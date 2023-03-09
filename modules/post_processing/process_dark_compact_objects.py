"""
This file contains post-processing to account for dim compact objects,
based on PopSyCLE (Lam et al. 2020).
"""
__all__ = ["ProcessDarkCompactObjects", ]
__author__ = "M.J. Houston"
__date__ = "2023-02-28"
__license__ = "GPLv3"
__version__ = "1.0.0"

import pandas
import numpy as np
from ._post_processing import PostProcessing


class ProcessDarkCompactObjects(PostProcessing):
    def __init__(self, model, **kwargs):
        """
        Parameters:
            model:
                SynthPop main model object
        """
        self.model = model
        self.logger = logger

    @staticmethod
    def mass_bh(m_zams, f_ej=0.9):
        """
        Black hole mass calculation with two branches split at initial mass of 42.21 solar masses
        Based on PopSyCLE (Lam et al., 2020) which draws from Raithel et al. 2018

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
        # branch i
        m_bh_core_i = -2.049 + 0.4140 * m_zams
        m_bh_all_i = 15.52 - 0.3294 * (m_zams - 25.97) - 0.02121 * (
                    m_zams - 25.97) ** 2 + 0.003120 * (m_zams - 25.97) ** 3
        # branch ii
        m_bh_core_ii = 5.697 + 7.8598 * 10 ** 8 * m_zams ** -4.858
        # branch determination: 0 for i and 1 for ii
        branch = (m_zams > 42.21).astype(int)
        # mass calculation
        return (f_ej * m_bh_core_i + (1 - f_ej) * m_bh_all_i) * (1 - branch) + m_bh_core_ii * branch

    @staticmethod
    def mass_ns(m_zams):
        """
        Neutron star final mass calculation, adopting the 1.6 solar mass average
        as a constant value for all initial masses
        Based on PopSyCLE (Lam et al., 2020) which draws from Raithel et al. 2018

        Parameters
        ----------
        m_zams
            float or array of float values for initial stellar mass in units of solar mass

        Returns
        -------
        m_ns
            float or array of float values for final neutron star mass in units of solar mass
        """
        return m_zams / m_zams * 1.6

    @staticmethod
    def mass_wd(m_zams):
        """
        White dwarf final mass calculation.
        Based on PopSyCLE (Lam et al., 2020) which draws from Raithel et al. 2018

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

    @staticmethod
    def compact_type(m_zams, phase):
        """
        Probabilistic drawing of compact object types
        Based on PopSyCLE (Lam et al., 2020) which draws from Raithel et al. 2018
        Probabilities from Lam et al. (2020) Table 1

        Parameters
        ----------
        m_zams
            array of float values for initial stellar mass in units of solar mass
        phase
            indicator that a star is in its final phase.
        Returns
        -------
        m_type
            array of integer values indicating object type
            0 = non-compact object or luminous white dwarf
            1 = dim white dwarf
            2 = neutron star
            3 = black hole
        """
        # Initial mass bin borders
        m_bins = np.array([0.5, 9, 15, 17.8, 18.5, 21.7, 25.2, 27.5, 60, 120])
        # probability for each type for each mass bin
        p_pre = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        p_wd = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        p_ns = np.array([0, 0, 1, 0.679, 0.833, 0.500, 0, 0.652, 0, 0.400])
        p_bh = np.array([0, 0, 0, 0.321, 0.167, 0.500, 1, 0.348, 1, 0.600])
        # arrange probabilities for drawing with numpy random
        ps = np.transpose([p_pre, p_wd, p_ns, p_bh])
        # probabilistic determination of compact object type for objects beyond
        return np.array(list(map(
            lambda m: np.random.choice(4, p=ps[np.searchsorted(m_bins, m)]), m_zams))) * phase

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        This is a placeholder for the postprocessing
        replace it with what ever you think is useful
        It must accept the pandas data frame and must return a pandas data frame

        Parameters
        ----------
        dataframe : dataframe
            original SynthPop output as pandas data frame

        Returns
        -------
        dataframe : dataframe
            modified pandas data frame
        """

        # determine object types for objects beyond isochrone grid
        m_type = self.compact_type(dataframe['iMass'], dataframe['In_Final_Phase'])
        # determine mass given object type
        m_compact = (self.mass_bh(dataframe['iMass']) * (m_type == 3).astype(int) +
                     self.mass_ns(dataframe['iMass']) * (m_type == 2).astype(int) +
                     self.mass_wd(dataframe['iMass']) * (m_type == 1).astype(int))
        # update final masses for objects beyond the isochrone grids
        m_final = (dataframe['In_Final_Phase'] * m_compact +
                   (1 - dataframe['In_Final_Phase']) * dataframe['Mass'])
        dataframe['Mass'] = m_final
        # add data column for object type
        dataframe['Dim_Compact_Object_Flag'] = m_type
        # get array indices for compact objects
        i_dim = m_type > 0
        if len(self.model.parms.opt_iso_props) > 0:
            i_dim = i_dim & np.isnan(dataframe[self.model.parms.col_names[0]])
        # set dim object magnitudes to nan
        for magcol in self.model.parms.chosen_bands:
            dataframe[magcol].loc[i_dim] = np.nan

        return dataframe
