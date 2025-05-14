"""
Post-processing to account for dim compact objects, based on PopSyCLE (Rose et al 2022).

The module will take all objects that have evolved past the MIST grid, and assign them a final
mass, removing their photometry. Optionally, one can just remove all of these objects instead.
"""

__all__ = ["ProcessDarkCompactObjects", ]
__author__ = "M.J. Huston"
__date__ = "2024-05-23"

import pandas
import numpy as np
from ._post_processing import PostProcessing

class ProcessDarkCompactObjects(PostProcessing):
    """
    Post-processing to account for dim compact objects, based on PopSyCLE (Rose et al 2022).

    Attributes
    ----------
    remove=False : boolean
        if true, remove dark compact objects from the catalog
    ifmr_name='SukhboldN20' : string
        selected initial-final mass relation;
        options are 'SukhboldN20', 'Raithel18', 'Spera15'
    """

    def __init__(self, model, logger, remove=False, ifmr_name='SukhboldN20', **kwargs):
        super().__init__(model, logger, **kwargs)
        self.remove = remove
        #: initial-final mass relation name to determine compact object masses.
        #: options are Raithel18, SukhboldN20, Spera15
        self.ifmr_name= ifmr_name

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
            m_bh_prelim = (1-branch)*m_bh_0 + branch*((1-f_z)*m_bh_0 + f_z*m_bh_zsun)
            # Assign any BHs < 3.0Msun to NS instead
            m_ns_backup = (m_bh_prelim<3.0).astype(int)*self.mass_ns(m_zams)
            m_bh = np.maximum(m_bh_prelim, m_ns_backup)
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
        z_cat_1 = np.array([(z>4.0e-3).astype(int), ((z<=4.0e-3)&(z>=1.0e-3)).astype(int), ((z<1.0e-3)).astype(int)])
        b1s = np.array([59.63    - 2.969e3*z  + 4.988e4*z**2,     40.98    + 3.415e4*z - 8.064e6*z**2,    np.repeat(67.07,    len(z))])
        k1s = np.array([45.04    - 2.176e3*z  + 3.806e4*z**2,     35.17    + 1.548e4*z - 3.759e6*z**2,    np.repeat(46.89,    len(z))])
        k2s = np.array([138.9    - 4.664e3*z  + 5.106e4*z**2,     20.36    + 1.162e5*z - 2.276e7*z**2,    np.repeat(113.8,    len(z))])
        d1s = np.array([2.790e-2 - 1.780e-2*z + 77.05*z**2,       2.500e-2 - 4.346*z   + 1.340e3*z**2,    np.repeat(2.199e-2, len(z))])
        d2s = np.array([6.730e-3 + 2.690*z    - 52.39*z**2,       1.750e-2 + 11.39*z   - 2.902e3*z**2,    np.repeat(2.602e-2, len(z))])

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
        if self.ifmr_name=='Raithel18':
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
        elif self.ifmr_name=='SukhboldN20':
            # Get value for metallicity dependence
            f_z = np.minimum(10**feh, np.ones(len(feh)))
            # Draw random numbers for bins that can be either NS or BH
            n_rand = np.random.uniform(size=len(m_zams))
            # Start with pre-CO objects, then go through mass bins, and assign appropriate type
            result = np.zeros(len(m_zams))
            result += ((m_zams>0.5)  & (m_zams<=9))    * 1
            result += ((m_zams>9)    & (m_zams<=15))   * 2
            result += ((m_zams>15)   & (m_zams<=21.8)) * ((n_rand<0.75)*2 + (n_rand>=0.75)*3)
            result += ((m_zams>21.8) & (m_zams<=25.2)) * 3
            result += ((m_zams>25.2) & (m_zams<=27.4)) * 2
            result += ((m_zams>27.4) & (m_zams<=60))   * 3
            result += ((m_zams>60)   & (m_zams<=120))  * ((n_rand<0.80*f_z)*2 + (n_rand>=0.80*f_z)*3)
        return result

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Perform the post-processing and return the modified DataFrame.
        """

        # Pick out which stars need processed
        in_final_phase = np.array(dataframe['In_Final_Phase'])
        proc_stars = dataframe[dataframe['In_Final_Phase']==1]
        # If we want to remove compact objects, do so and return
        if self.remove:
            return dataframe.drop(proc_stars.index)
        # Otherwise, we need to handle the compact objects properly. 
        m_init = np.array(dataframe['iMass'])
        feh_init = np.array(dataframe['Fe/H_initial'])
        m_final_pre = np.array(dataframe['Mass'])
        
        # For IFMRs with probabilistic object types
        if self.ifmr_name in ['Raithel18', 'SukhboldN20']:
            # Probabilistic determination of object types
            m_type = self.compact_type_from_initial(m_init, feh_init)
            # Get possible masses and select by type
            m_compact = (self.mass_bh(m_init, feh_init) * (m_type == 3).astype(int) +
                         self.mass_ns(m_init) * (m_type == 2).astype(int) +
                         self.mass_wd(m_init) * (m_type == 1).astype(int))
            m_type = self.compact_type_from_final(m_compact)
        # For IFMRs with analytic mass determination, then type assigned by mass
        elif self.ifmr_name in ['Spera15']:
            # Cycle through evolved stars, calculating mass & type
            m_compact = self.mass_spera15(m_init,feh_init)
            m_type = self.compact_type_from_final(m_compact)
        
        # Results into data frame
        m_final = (in_final_phase * m_compact +
                   (1 - in_final_phase) * m_final_pre)
        dataframe['Mass'] = m_final
        dataframe['Dim_Compact_Object_Flag'] = m_type*in_final_phase

        # Set dim object magnitudes to nan
        for magcol in self.model.parms.chosen_bands:
            dataframe.loc[proc_stars.index, magcol] = np.nan
            
        return dataframe
