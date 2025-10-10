"""
Evolution module to store information for the SpiseaGenerator. Not valid for a standard StarGenerator.
"""

__all__ = ["SpiseaCluster", ]
__author__ = "M.J. Huston"
__date__ = "2025-05-28"

from ._evolution import EvolutionIsochrones, EvolutionInterpolator, EVOLUTION_DIR
import numpy as np
import json

class SpiseaCluster(EvolutionIsochrones,EvolutionInterpolator):
    """
    Placeholder object to store settings for use by the SpiseaGenerator, which 
    will generate and evolve stars as binned SPISEA clusters.
    """
    def __init__(self, columns, evo_model_name="MISTv1.2-synthpop", atm_func_name="get_merged_atmosphere", 
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
        with open(f"{EVOLUTION_DIR}/spisea_effective_wavelengths.json") as f:
            self.eff_wavelengths = json.load(f)


    def get_cols(self, columns):
        with open(f"{EVOLUTION_DIR}/spisea_filters.json") as f:
            all_spisea_filters = json.load(f)
        magsys = {}
        all_bands = []
        non_mag_cols = []

        for column in columns:
            column_split = column.split('-')
            # Non-magnitude columns
            if column in self.allowed_non_mag_cols:
                non_mag_cols.append(column)
            # Magnitude systems
            elif column in all_spisea_filters.keys():
                magsys = all_spisea_filters[column]
                # Magnitude systems with nested categories & need all
                if hasattr(magsys, "keys"):
                    for magsys_subset in magsys.keys():
                        for band in magsys[magsys_subset]:
                            all_bands.append(column+'-'+magsys_subset+'-'+band)
                # Magnitude systems with no nested catagories & need all
                else:
                    for band in magsys:
                        all_bands.append(column+'-'+band)
            # Magnitude systems with nested or band selections
            elif (len(column_split) in [2,3]) and (column_split[0] in all_spisea_filters.keys()):
                magsys = all_spisea_filters[column_split[0]]
                # Specified nested band set, use all filters
                if (len(column_split)==2) and hasattr(magsys, "keys"):
                    if column_split[1] in magsys.keys():
                        for band in magsys[column_split[1]]:
                            all_bands.append(column+'-'+band)
                    else:
                        raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                # Specified nested band set, select filters
                elif len(column_split)==3:
                    if column_split[1] in magsys.keys():
                        if column_split[2] in magsys[column_split[1]]:
                            all_bands.append(column)
                        else:
                            raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                    else:
                        raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                elif len(column_split)==2:
                    if column_split[1] in magsys:
                        all_bands.append(column)
                    else:
                        raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
                else:
                    raise ValueError('Invalid column '+column+' for SPISEA isochrones.')
            else:
                raise ValueError('Invalid column '+column+' for SPISEA isochrones.')

        return None, list(np.unique(non_mag_cols)), list(np.unique(all_bands))

    def get_evolved_props(**kwargs):
        return
