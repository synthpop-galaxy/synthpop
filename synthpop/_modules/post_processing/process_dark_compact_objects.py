"""
This file contains post-processing to account for dim compact objects,
based on PopSyCLE (Rose et al 2022).
"""
__all__ = ["ProcessDarkCompactObjects", ]
__author__ = "M.J. Huston"
__date__ = "2024-05-23"
__license__ = "GPLv3"
__version__ = "1.0.0"

import pandas
import numpy as np
from ._post_processing import PostProcessing

#TODO finish implemenmtation of SukhboldN20 and Spera15 options
# Currently, SukhboldN20 works with a simplification assuming solar metallicity 
# when determining compact object type
# Spera15 is theoretically working
# All need tested after completion
class ProcessDarkCompactObjects(PostProcessing):
    def __init__(self, model, logger, remove=False, ifmr_name='Raithel18', **kwargs):
        """
        Parameters:
            model:
                SynthPop main model object
    	    remove:
    		  if True, remove compact objects from catalog; if False, keep them in
    	    ifmr_name:
    		  select initial final mass relation to determine compact object masses based on their initial mass
    		  current options are Raithel18, SukhboldN20
        """
        self.model = model
        self.remove = remove
        self.ifmr_name= ifmr_name
        self.logger = logger

    def mass_bh(self, m_zams, feh, f_ej=0.9):
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
        if self.ifmr_name=='Raithel18':
            m_bh_core_i = -2.049 + 0.4140 * m_zams
            m_bh_all_i = 15.52 - 0.3294 * (m_zams - 25.97) - 0.02121 * (
                        m_zams - 25.97) ** 2 + 0.003120 * (m_zams - 25.97) ** 3
            # branch ii
            m_bh_core_ii = 5.697 + 7.8598 * 10 ** 8 * m_zams ** -4.858
            # branch determination: 0 for i and 1 for ii
            branch = (m_zams > 42.21).astype(int)
            m_bh = (f_ej * m_bh_core_i + (1 - f_ej) * m_bh_all_i) * (1 - branch) + m_bh_core_ii * branch
        elif self.ifmr_name=='SukhboldN20':
            f_z = np.minimum(10**feh, np.ones(len(feh)))
            m_bh_0 = 0.4652*m_zams - 3.2917
            m_bh_zsun = -0.271*m_zams + 24.743
            branch = (m_zams > 39.6).astype(int)
            m_bh = (1-branch)*m_bh_0 + branch*((1-f_z)*m_bh_0 + f_z*m_bh_zsun)
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

    def single_mass_spera15(self, m_zams, feh):
        """
        Remnant mass calculation from Spera et al 2015, appendix C
        Takes in the m_zams and Fe/H for a single star at a time
        """
        # First, calculate M_CO, based on M_ZAMS and Z
        # Note: z equation from Rose et al 2022, 
        #       with z_sun from Ekstroem et al 2012
        z = 0.014*10**feh
        if z>4.0e-3:
            # C13
            b1 = 59.63 - 2.969e3*z + 4.988e4*z**2
            k1 = 45.04 - 2.176e3*z + 3.806e4*z**2
            k2 = 138.9 - 4.664e3*z + 5.106e4*z**2
            d1 = 2.790e-2 - 1.780e-2*z + 77.05*z**2
            d2 = 6.730e-3 + 2.690*z - 52.39*z**2
        elif z>1.0e-3:
            # C14
            b1 = 40.98 + 3.415e4*z - 8.064e6*z**2
            k1 = 35.17 + 1.548e4*z - 3.759e6*z**2
            k2 = 20.36 + 1.162e5*z - 2.276e7*z**2
            d1 = 2.500e-2 - 4.346*z + 1.340e3*z**2
            d2 = 1.750e-2 + 11.39*z - 2.902e3*z**2
        else:
            # C15
            b1 = 67.07
            k1 = 46.89
            k2 = 113.8
            d1 = 2.199e-2
            d2 = 2.602e-2
        # C12
        g1 = 0.5 / (1 + 10**((k1-m_zams)*d1))
        g2 = 0.5 / (1 + 10**((k2-m_zams)*d2))
        # C11
        m_co = -2.0 + (b1+2.0)*(g1+g2)
        # Then, we calculate M_REM
        if z<= 5e-4:
            # C3
            m = -6.476e2*z + 1.911
            q = 2.300e3*z + 11.67
            # C2
            p = -2.333 + 0.1559*m_co + 0.2700*m_co**2
            f = m*m_co + q
            # C1
            if m_co < 5:
                m_rem = max(p, 1.27)
            elif m_co<10:
                m_rem = p 
            else:
                m_rem = min(p, f)
        else:
            if z>2.0e-3:
                # C8
                m = 1.217 
                q = 1.061
            elif z>1.0e-3:
                # C9
                m = -43.82*z + 1.340
                q = -1.296e4*z + 26.98 
            else:
                # C10
                m = -6.476e2*z + 1.911 
                q = 2.300e3*z + 11.67
            if z>1e-3:
                # C6
                a1 = 1.340 - 29.46 / (1 + (z/1.110e-3)**2.361)
                a2 = 80.22 - 74.73 * z**0.965 / (2.720e-3 + z**0.965)
                l = 5.683 + 3.533 / (1 + (z/7.430e-3)**1.993)
                eta = 1.066 - 1.121 / (1 + (z/2.558e-2)**0.609)
            else:
                # C7
                a1 = 1.105e5*z - 1.258e2
                a2 = 91.56 - 1.957e4*z - 1.558e7*z**2
                l = 1.134e4*z - 2.143 
                eta = 3.090e-2 - 22.30*z + 7.363e4*z**2
            # C5
            h = a1 + (a2-a1)/(1+10**((l-m_co)*eta))
            f = m*m_co+q
            # C4
            if m_co<5:
                m_rem = max(h, 1.27)
            elif m_co<10:
                m_rem = h 
            else:
                m_rem = max(h, f)
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
            nteger value indicating object type
            1 = dim white dwarf
            2 = neutron star
            3 = black hole
        """
        if m_fin<1.4:
            return 1
        elif m_fin<3:
            return 2
        else:
            return 3

    def compact_type_from_initial(self, m_zams, phase):
        """
        Probabilistic drawing of compact object types
        Based on PopSyCLE (Lam et al 2020; Rose et al 2022)
        Which pulls from Rathiel et al 2018 and Sukhbold et al 2020

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
        if self.ifmr_name=='Raithel18':
            # Initial mass bin borders
            m_bins = np.array([0.5, 9, 15, 17.8, 18.5, 21.7, 25.2, 27.5, 60, 120])
            # probability for each type for each mass bin
            p_pre = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            p_wd = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
            p_ns = np.array([0, 0, 1, 0.679, 0.833, 0.500, 0, 0.652, 0, 0.400])
            p_bh = np.array([0, 0, 0, 0.321, 0.167, 0.500, 1, 0.348, 1, 0.600])
            # arrange + draw from probabilities
            ps = np.transpose([p_pre, p_wd, p_ns, p_bh])
            result = np.array(list(map(
                lambda m: np.random.choice(4, p=ps[np.searchsorted(m_bins, m)]), m_zams))) * phase
        elif self.ifmr_name=='SukhboldN20':
            # Initial mass bin borders
            m_bins = np.array([0.5, 9, 15, 21.8, 25.2, 27.4, 60, 120])
            # probability for each type for each mass bin
            p_pre = np.array([1, 0, 0, 0, 0, 0, 0, 0])
            p_wd = np.array([0, 1, 0, 0, 0, 0, 0, 0])
            p_ns = np.array([0, 0, 1, 0.75, 0, 1, 0, 0.8])
            p_bh = np.array([0, 0, 0, 0.25, 1, 0, 1, 0.2]) #TODO this is cheating the metallicity rn - need to fix it
            # arrange + draw from probabilities
            ps = np.transpose([p_pre, p_wd, p_ns, p_bh])
            result = np.array(list(map(
                lambda m: np.random.choice(4, p=ps[np.searchsorted(m_bins, m)]), m_zams))) * phase
        # probabilistic determination of compact object type for objects beyond
        return result 

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Function that executes the post-processing to assign compact object masses

        Parameters
        ----------
        dataframe : dataframe
            original SynthPop output as pandas data frame

        Returns
        -------
        dataframe : dataframe
            modified pandas data frame
        """

        # Pick out which stars need processed
        proc_stars = dataframe[dataframe['In_Final_Phase']==1]
        # If we want to remove compact objects, do so and return
        if self.remove:
            return dataframe.drop(proc_stars.index)
        # Otherwise, we need to handle the compact objects properly. 
        # Start by adding a column to flag the compact objects
        dataframe['Dim_Compact_Object_Flag'] = np.zeros(len(dataframe.index))

        # For IFMRs with probabilistic object types
        if self.ifmr_name in ['Raithel18', 'SukhboldN20']:
            # Probabilistic determination of object types
            m_type = self.compact_type_from_initial(dataframe['iMass'], dataframe['In_Final_Phase'])
            # Get possible masses and select by type
            m_compact = (self.mass_bh(dataframe['iMass'], dataframe['Fe/H_initial']) * (m_type == 3).astype(int) +
                         self.mass_ns(dataframe['iMass']) * (m_type == 2).astype(int) +
                         self.mass_wd(dataframe['iMass']) * (m_type == 1).astype(int))
            # Update masses in data frame
            m_final = (dataframe['In_Final_Phase'] * m_compact +
                       (1 - dataframe['In_Final_Phase']) * dataframe['Mass'])
            dataframe['Mass'] = m_final
            # Add flags for compact objects
            dataframe['Dim_Compact_Object_Flag'] = m_type
        # For IFMRs with analytic mass determination, then type assigned by mass
        elif self.ifmr_name in ['Spera15']:
            # Cycle through evolved stars, calculating mass & type
            for idx in proc_stars.index:
                m_final = self.single_mass_spera15(proc_stars['iMass'][idx],
                    proc_stars['Fe/H_initial'][idx])
                m_type = self.compact_type_from_final(m_final)
                dataframe['Mass'][idx] = m_final  
                dataframe['Dim_Compact_Object_Flag'][idx] = m_type

        # Set dim object magnitudes to nan
        for magcol in self.model.parms.chosen_bands:
            dataframe[magcol].loc[proc_stars.index] = np.nan
            
        return dataframe
