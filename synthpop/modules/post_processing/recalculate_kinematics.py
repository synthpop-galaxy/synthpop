"""
Post-processing module to recalculate and replace the kinematics of a catalog.
"""

__all__ = ["RecalculateKinematics", ]
__author__ = "M.J. Huston"
__date__ = "2025-05-14"

import pandas
import numpy as np
from ._post_processing import PostProcessing

class RecalculateKinematics(PostProcessing):
    """
    Post-processing module to recalculate and replace the kinematics of a catalog.
    
    Attributes
    ----------
    pop_ids : list
        Selection of which populations to replace kinematics for.
        If none is given, all will be replaced.
    """

    def __init__(self, model, logger, pop_ids=None, **kwargs):
        super().__init__(model,logger, **kwargs)
        self.pop_ids = pop_ids

    def do_post_processing(self, dataframe: pandas.DataFrame) -> pandas.DataFrame:
        """
        Replace the kinematic columns and return the modified catalog.
        """
        if self.pop_ids is not None:
            pop_ids = self.pop_ids
        else:
            pop_ids = np.unique(dataframe['pop'].to_numpy())


        for pop in pop_ids:
            idxs = dataframe[dataframe['pop']==float(pop)].index.to_numpy()
            u, v, w, vr, mu_l, mu_b, vr_lsr = self.model.populations[int(pop)].do_kinematics(
                                                dataframe.Dist[idxs].to_numpy(), dataframe.l[idxs].to_numpy(), dataframe.b[idxs].to_numpy(),                    
                                                dataframe.x[idxs].to_numpy(), dataframe.y[idxs].to_numpy(), dataframe.x[idxs].to_numpy())
            dataframe.loc[idxs, 'U'] = u
            dataframe.loc[idxs, 'V'] = v
            dataframe.loc[idxs, 'W'] = w
            dataframe.loc[idxs, 'mul'] = mu_l
            dataframe.loc[idxs, 'mub'] = mu_b
            dataframe.loc[idxs, 'vr_bc'] = vr
            dataframe.loc[idxs, 'VR_LSR'] = vr_lsr

        return dataframe
