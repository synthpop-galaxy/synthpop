"""
Post-processing to add random kick velocities to neutron stars and black holes, according to 
a Maxwellian distribution with a user-input mean.
"""

__all__ = ["NatalKicks", ]
__author__ = "M.J. Huston"
__date__ = "2025-10-15"

import pandas as pd
import numpy as np
from ._post_processing import PostProcessing
from scipy.stats import maxwell, uniform_direction
from synthpop.synthpop_utils.coordinates_transformation import CoordTrans
import pdb

class NatalKicks(PostProcessing):
    """
    Post-processing to add kicks to NSs and BHs, based on PopSyCLE (Rose et al 2022).

    Attributes
    ----------
    kick_mean_bh=0 : float
        mean of the maxwellian kick distribution for black holes (km/s)
    kick_mean_ns=0 : float
        mean of the maxwellian kick distribution for neutron stars (km/s)
    """

    def __init__(self, model, logger, kick_mean_bh=0, kick_mean_ns=0, **kwargs):
        super().__init__(model, logger, **kwargs)
        self.kick_mean_ns = kick_mean_ns
        self.kick_mean_bh = kick_mean_bh

    def do_post_processing(self, system_df: pd.DataFrame,
            companion_df: pd.DataFrame):
        """
        Perform the post-processing and return the modified DataFrame.
        """

        # Pick out which stars need processed
        phase = system_df['phase'].to_numpy().astype(int)
        proc_stars = system_df.index

        # Apply birth kicks
        kick_idxs = proc_stars[phase>=102]
        kick_mtypes = phase[phase>=102]
        # Generate random velocities
        kick_vel = maxwell.rvs(size=len(kick_idxs), scale=1, loc=0) * \
                            self.kick_mean_ns**(kick_mtypes==102).astype(int) * \
                            self.kick_mean_bh**(kick_mtypes==103).astype(int)
        # Generate random directions
        rand_dir = uniform_direction.rvs(dim=3, size=len(kick_idxs))
        # Update cartesian velocities
        l_deg = system_df['l'][kick_idxs].to_numpy()
        b_deg = system_df['b'][kick_idxs].to_numpy()
        u_new = system_df['U'][kick_idxs].to_numpy() + kick_vel * rand_dir[:,0]
        v_new = system_df['V'][kick_idxs].to_numpy()  + kick_vel * rand_dir[:,1]
        w_new = system_df['W'][kick_idxs].to_numpy()  + kick_vel * rand_dir[:,2]
        coord_trans = CoordTrans(sun=self.model.parms.sun)
        system_df.loc[kick_idxs, 'U'] = u_new
        system_df.loc[kick_idxs, 'V'] = v_new
        system_df.loc[kick_idxs, 'W'] = w_new
        # Convert to and update proper motion/radial velocities
        kick_ls = system_df['l'][kick_idxs].to_numpy()
        kick_bs = system_df['b'][kick_idxs].to_numpy()
        kick_dists = system_df['Dist'][kick_idxs].to_numpy()
        if 'mul' in system_df:
            vr_new, mul_new, mub_new = coord_trans.uvw_to_vrmulb(kick_ls, kick_bs, kick_dists, u_new, v_new, w_new)
            system_df.loc[kick_idxs, 'vr_bc'] = vr_new
            system_df.loc[kick_idxs, 'mul'] = mul_new
            system_df.loc[kick_idxs, 'mub'] = mub_new
        if 'mura' in system_df:
            vr_new, mura_new, mudec_new = coord_trans.uvw_to_vrmulb(kick_ls, kick_bs, kick_dists, u_new, v_new, w_new)
            system_df.loc[kick_idxs, 'vr_bc'] = vr_new
            system_df.loc[kick_idxs, 'mura'] = mura_new
            system_df.loc[kick_idxs, 'mudec'] = mudec_new
        system_df.loc[kick_idxs, 'VR_LSR'] = coord_trans.vr_bc_to_vr_lsr(l_deg,b_deg,vr_new)
            
        return system_df, companion_df
