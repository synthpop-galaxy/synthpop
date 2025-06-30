"""
Postprocessing module to convert from galactic coordinates to equatorial.
"""

__all__ = ["EquatorialCoordinates", ]
__author__ = "M.J. Huston"
__date__ = "2025-06-30"

import pandas
import numpy as np
from ._post_processing import PostProcessing

class EquatorialCoordinates(PostProcessing):
    """
    Postprocessing module to convert magnitude systems for any
    filters provided by MIST. Allowed systems are Vega, AB, and ST.
    
    Attributes
    ----------
    keep_galactic : bool
        selection for whether to keep galactic coordinates in data frame or remove
    """

    def __init__(self, model, logger, keep_galactic=False, **kwargs):
        super().__init__(model,logger, **kwargs)
        self.keep_galactic = keep_galactic

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Perform the magnitude conversions and returns the modified DataFrame.
        """
        # Create new columns where needed
        l_idx = np.where(dataframe.columns=='l')[0][0]
        dataframe.insert(l_idx, 'ra', np.nan)
        dataframe.insert(l_idx+1, 'dec', np.nan)
        try:
            mul_idx = np.where(dataframe.columns=='mul')[0][0]
            dataframe.insert(mul_idx, 'mura', np.nan)
            dataframe.insert(mul_idx+1, 'mudec', np.nan)
            no_kinematics=False
        except:
            no_kinematics = True

        # Loop over populations to convert coordinate systems and fill in dataframe column
        for popid in range(len(self.model.populations)):
            stars_in_pop = np.where(dataframe['pop']==float(popid))[0]
            l_pop, b_pop = dataframe['l'][stars_in_pop].to_numpy(), dataframe['b'][stars_in_pop].to_numpy()
            dist_pop = dataframe['Dist'][stars_in_pop].to_numpy()
            ra_pop, dec_pop = self.model.populations[popid].coord_trans.lb_to_ad(l_pop,b_pop)
            dataframe.loc[stars_in_pop, 'ra'] = ra_pop
            dataframe.loc[stars_in_pop, 'dec'] = dec_pop
            if not no_kinematics:
                u_pop, v_pop, w_pop = dataframe['U'][stars_in_pop].to_numpy(), dataframe['V'][stars_in_pop].to_numpy(), \
                                        dataframe['W'][stars_in_pop].to_numpy()
                vr_pop, mura_pop, mudec_pop = self.model.populations[popid].coord_trans.uvw_to_vrmuad(l_pop,b_pop,dist_pop, u_pop,v_pop,w_pop)
                dataframe.loc[stars_in_pop, 'mura'] = mura_pop
                dataframe.loc[stars_in_pop, 'mudec'] = mudec_pop

        # Remove galactic coordinates if selected
        if not self.keep_galactic:
            dataframe.drop(columns=['l','b'],inplace=True)
            if not no_kinematics:
                dataframe.drop(columns=['mul','mub'],inplace=True)
        
        return dataframe
