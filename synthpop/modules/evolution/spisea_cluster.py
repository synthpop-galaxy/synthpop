"""
Evolution module to store information for the SpiseaGenerator. Not valid for a standard StarGenerator.
"""

__all__ = ["SpiseaCluster", ]
__author__ = "M.J. Huston"
__date__ = "2025-05-28"

from ._evolution import EvolutionIsochrones, EvolutionInterpolator, EVOLUTION_DIR
import numpy as np
import json
from spisea import evolution as spisea_evolution
from spisea import atmospheres as spisea_atmospheres
from spisea import synthetic as spisea_synthetic

class SpiseaCluster(EvolutionIsochrones,EvolutionInterpolator):
    """
    Placeholder object to store modules for use by the SpiseaGenerator, which 
    will generate and evolve stars as binned SPISEA clusters.
    """
    def __init__(self, columns, spisea_evolution_name="MISTv1", block_spisea_prints=True,
                    spisea_evolution_kwargs={"version":1.2, "synthpop_extension":True},
                    spisea_atm_func_name="get_merged_atmosphere", spisea_wd_atm_func_name="get_wd_atmosphere",
                    multiplicity_name=None, min_mass=0, max_mass=1000, **kwargs):
        self.block_spisea_prints=block_spisea_prints
        if spisea_evolution_name=='MISTv1':
            self.feh_list = np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.25,0.5])
            self.log_age_list = np.linspace(5.0,10.3,107)
            self.log_age_list[0] = 5.01
        elif spisea_evolution_name=='MergedBaraffePisaEkstromParsec':
            self.feh_list = np.array([0.0])
            self.log_age_list = np.linspace(6.0,10.1,83)
            self.log_age_list[-1] = 10.9
            Warning(f"Evolution module {spisea_evolution_name} only includes solar metallicy. All stars will be assigned solar metallicity.")
        else:
            raise ValueError("Invalid SPISEA evolution_model. Only MISTv1 and MergedBaraffePisaEkstromParsec are available at this time.")
        self.spisea_evolution = getattr(spisea_evolution, spisea_evolution_name)(**spisea_evolution_kwargs)

        self.spisea_atm_func = getattr(spisea_atmospheres, spisea_atm_func_name)
        self.spisea_wd_atm_func = getattr(spisea_atmospheres, spisea_wd_atm_func_name)
        self.min_mass = min_mass
        self.max_mass = max_mass

        self.allowed_non_mag_cols = ["[Fe/H]_init", "log10_isochrone_age_yr", 'phase',
            'star_mass', 'initial_mass', 'log_L', 'log_R', 'log_Teff', 'log_g']

        self.magsys, self.non_mag_cols, self.bands, self.bands_obs_str = self.get_cols(columns)
        with open(f"{EVOLUTION_DIR}/spisea_effective_wavelengths.json") as f:
            all_eff_wavelengths = json.load(f)
        self.eff_wavelengths = {self.bands[i]:all_eff_wavelengths[self.bands_obs_str[i]] for i in range(len(self.bands))}


    def get_cols(self, columns):
        with open(f"{EVOLUTION_DIR}/spisea_filters.json") as f:
            all_spisea_filters = json.load(f)
        magsys = {}
        all_bands = []
        non_mag_cols = []

        for column in columns:
            column_split = column.split(',')
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
                            all_bands.append(column+','+magsys_subset+','+band)
                # Magnitude systems with no nested catagories & need all
                else:
                    for band in magsys:
                        all_bands.append(column+','+band)
            # Magnitude systems with nested or band selections
            elif (len(column_split) in [2,3]) and (column_split[0] in all_spisea_filters.keys()):
                magsys = all_spisea_filters[column_split[0]]
                # Specified nested band set, use all filters
                if (len(column_split)==2) and hasattr(magsys, "keys"):
                    if column_split[1] in magsys.keys():
                        for band in magsys[column_split[1]]:
                            all_bands.append(column+','+band)
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

        all_bands_obs_str = list(np.unique(all_bands))
        all_bands = ['m_'+spisea_synthetic.get_filter_col_name(band) for band in all_bands_obs_str]

        return None, list(np.unique(non_mag_cols)), all_bands, all_bands_obs_str

    def get_evolved_props(**kwargs):
        raise ValueError("SpiseaCluster Evolution module is only compatible with SpiseaGenerator, not StarGenerator")
        return
