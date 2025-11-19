"""
Post-processing module to recalculate and replace the kinematics of a catalog.
"""

__all__ = ["RecalculateKinematics", ]
__author__ = "M.J. Huston"
__date__ = "2025-05-14"

import pandas as pd
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

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Replace the kinematic columns and return the modified catalog.
        """
        if self.pop_ids is not None:
            pop_ids = self.pop_ids
        else:
            pop_ids = np.unique(system_df['pop'].to_numpy())


        for pop in pop_ids:
            idxs = system_df[system_df['pop']==float(pop)].index.to_numpy()
            u, v, w, vr, mu_l, mu_b, vr_lsr = self.model.populations[int(pop)].do_kinematics(
                                                system_df.Dist[idxs].to_numpy(), system_df.l[idxs].to_numpy(), system_df.b[idxs].to_numpy(),
                                                system_df.x[idxs].to_numpy(), system_df.y[idxs].to_numpy(), system_df.x[idxs].to_numpy())
            system_df.loc[idxs, 'U'] = u
            system_df.loc[idxs, 'V'] = v
            system_df.loc[idxs, 'W'] = w
            system_df.loc[idxs, 'mul'] = mu_l
            system_df.loc[idxs, 'mub'] = mu_b
            system_df.loc[idxs, 'vr_bc'] = vr
            system_df.loc[idxs, 'VR_LSR'] = vr_lsr

        return system_df, companion_df
