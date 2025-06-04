"""
Evolution module to store information for the SpiseaGenerator. Not valid for a standard StarGenerator.
"""

__all__ = ["Spisea", ]
__author__ = "M.J. Huston"
__date__ = "2024-11-20"
__license__ = "GPLv3"
__version__ = "1.0.0"

from ._evolution import EvolutionIsochrones, EvolutionInterpolator
import numpy as np

class SpiseaCluster(EvolutionIsochrones,EvolutionInterpolator):
    """
    Placeholder object to store settings for use by the SpiseaGenerator, which 
    will generate and evolve stars as binned SPISEA clusters.
    """
    def __init__(self, columns, evo_model_name="MISTv1.2", atm_func_name="get_merged_atmosphere", 
                    wd_atm_func_name="get_wd_atmosphere", ifmr_name="IFMR_N20_Sukhbold", 
                    multiplicity_name=None, min_mass=0, max_mass=1000, **kwargs):
        self.evo_model_name = evo_model_name
        self.atm_func_name = atm_func_name
        self.wd_atm_func_name = wd_atm_func_name
        self.ifmr_name=ifmr_name
        self.multiplicity_name = multiplicity_name
        self.min_mass = min_mass
        self.max_mass = max_mass

        self.allowed_non_mag_cols = ["[Fe/H]_init", "log10_isochrone_age_yr", 'phase',
            'star_mass', 'initial_mass', 'log_L', 'log_R', 'log_Teff', 'log_g']

        self.magsys, self.non_mag_cols, self.bands = self.get_cols(columns)


    def get_cols(self, columns):
        magsys = {}
        all_bands = []
        non_mag_cols = []

        for column in columns:
            if column in self.allowed_non_mag_cols:
                non_mag_cols.append(column)
            else:
                try:
                    column_split = column.split('-')
                    if len(column_split)==2:
                        msys = column_split[0]
                    elif len(column_split)==3:
                        msys = column_split[0]+'-'+column_split[1]
                    filt = column_split[-1]
                    if (msys not in magsys):
                        magsys[msys] = []
                    magsys[msys].append(column)
                    all_bands.append(column)
                except:
                    raise ValueError('Invalid column '+column+' for SPISEA isochrones.')

        return magsys, list(np.unique(non_mag_cols)), all_bands

    def get_evolved_props(**kwargs):
        return