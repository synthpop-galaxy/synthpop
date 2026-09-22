"""
Post-processing to add random kick velocities to neutron stars and black holes.
"""

__all__ = ["NatalKicks", ]
__author__ = "M.J. Huston"
__date__ = "2025-10-15"

import pandas as pd
import numpy as np
from ._post_processing import PostProcessing
import scipy.stats
from synthpop.synthpop_utils.coordinates_transformation import CoordTrans
import pdb

class NatalKicks(PostProcessing):
    """
    Post-processing to add kicks to NSs and BHs from any scipy.stats distribution
    with defaults from PopSyCLE.
    
    Note: scipy.stats distribution inputs can be non-intiutive, so be sure to check that
    your input quantities are scaled correctly.

    Attributes
    ----------
    ns_distribution : str
        name of scipy.stats random variable class for neutron stars
    ns_kwargs : dict
        arguments and keyword arguments by name for the indicated distribution
    bh_distribution : str
        name of scipy.stats random variable class for black holes
    bh_kwargs : dict
        arguments and keyword arguments by name for the indicated distribution
    """

    def __init__(self, model, logger, ns_distribution='lognorm',
                 ns_kwargs={'loc':0, 'scale': np.exp(5.6), 's': 0.68},
                 bh_distribution='maxwell',
                 bh_kwargs={'loc':0, 'scale': 100/(2*np.sqrt(2/np.pi))},
                 **kwargs):
        super().__init__(model, logger, **kwargs)
        self.ns_distribution = ns_distribution
        self.ns_kwargs = ns_kwargs
        self.bh_distribution = bh_distribution
        self.bh_kwargs = bh_kwargs

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Perform the post-processing and return the modified DataFrame.
        """
        if 'kick_x' in system_df:
            self.logger.log("Natal kicks already handled. Skipping postproc kicks.")
            return system_df, companion_df

        # Pick out which stars need processed
        phase = system_df['phase'].to_numpy().astype(int)
        proc_stars = system_df.index
        kick_idxs = proc_stars[phase>=102]

        # Get NS kicks
        ns_idxs = proc_stars[phase==102]
        ns_distr = getattr(scipy.stats, self.ns_distribution)
        ns_kick_vel = ns_distr.rvs(size=len(ns_idxs), **self.ns_kwargs)
        rand_dir = scipy.stats.uniform_direction.rvs(dim=3, size=len(ns_idxs))
        system_df.loc[ns_idxs, 'U'] += ns_kick_vel * rand_dir[:,0]
        system_df.loc[ns_idxs, 'V'] += ns_kick_vel * rand_dir[:,1]
        system_df.loc[ns_idxs, 'W'] += ns_kick_vel * rand_dir[:,2]
        system_df.loc[ns_idxs, 'kick_vel'] = ns_kick_vel

        # Get BH kicks
        bh_idxs = proc_stars[phase==103]
        bh_distr = getattr(scipy.stats, self.bh_distribution)
        bh_kick_vel = bh_distr.rvs(size=len(bh_idxs), **self.bh_kwargs)
        rand_dir = scipy.stats.uniform_direction.rvs(dim=3, size=len(bh_idxs))
        system_df.loc[bh_idxs, 'U'] += bh_kick_vel * rand_dir[:,0]
        system_df.loc[bh_idxs, 'V'] += bh_kick_vel * rand_dir[:,1]
        system_df.loc[bh_idxs, 'W'] += bh_kick_vel * rand_dir[:,2]
        system_df.loc[bh_idxs, 'kick_vel'] = bh_kick_vel

        # Convert to and update proper motion/radial velocities
        u_new = system_df['U'][kick_idxs].to_numpy()
        v_new = system_df['V'][kick_idxs].to_numpy()
        w_new = system_df['W'][kick_idxs].to_numpy()
        kick_ls = system_df['l'][kick_idxs].to_numpy()
        kick_bs = system_df['b'][kick_idxs].to_numpy()
        kick_dists = system_df['Dist'][kick_idxs].to_numpy()
        coord_trans = CoordTrans(sun=self.model.parms.sun)
        vr_new, mul_new, mub_new = coord_trans.uvw_to_vrmulb(kick_ls, kick_bs, kick_dists,
                                                             u_new, v_new, w_new)
        system_df.loc[kick_idxs, 'vr_bc'] = vr_new
        system_df.loc[kick_idxs, 'mul'] = mul_new
        system_df.loc[kick_idxs, 'mub'] = mub_new
        system_df.loc[kick_idxs, 'VR_LSR'] = coord_trans.vr_bc_to_vr_lsr(kick_ls, kick_bs, vr_new)
            
        return system_df, companion_df
