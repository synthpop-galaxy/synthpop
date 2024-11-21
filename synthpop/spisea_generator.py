
"""
This file contains the in-progress Spisea alternative to StarGenerator,
It bins stars by age and metallicity and generates SPISEA clusters,
    then assigns them locations based on the population_density
"""

__all__ = ["SpiseaGenerator"]
__author__ = "M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2024-11-20"
__license__ = "GPLv3"
__version__ = "0.2.0"

from typing import Set, Tuple, Dict
import numpy as np
import time
import spisea
import os

# Local Imports
# used to allow running as main and importing to another script
try:
    from . import constants as const
except ImportError:
    import constants as const
    import synthpop_utils as sp_utils
    from position import Position
    from synthpop_utils.synthpop_logging import logger
    from synthpop_utils import coordinates_transformation as coord_trans
    from synthpop_utils import Parameters
    from modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from modules.age import Age
    from modules.initial_mass_function import InitialMassFunction
    from modules.kinematics import Kinematics
    from modules.metallicity import Metallicity
    from modules.population_density import PopulationDensity
    from star_generator import StarGenerator

else:  # continue import when if synthpop is imported
    from . import synthpop_utils as sp_utils
    from .star_generator import StarGenerator
    from .position import Position
    from .synthpop_utils.synthpop_logging import logger
    from .synthpop_utils import coordinates_transformation as coord_trans
    from .synthpop_utils import Parameters
    from .modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from .modules.age import Age
    from .modules.initial_mass_function import InitialMassFunction
    from .modules.kinematics import Kinematics
    from .modules.metallicity import Metallicity
    from .modules.population_density import PopulationDensity

class SpiseaGenerator(StarGenerator):
    def __init__(self, imf_module, age_module, met_module, evolution_module,
            glbl_params, position, max_mass, logger, spisea_dir=const.DATA_DIR+'/isochrones/spisea/'):
        # General synthpop things
        self.age_module = age_module
        self.met_module = met_module
        self.kinematics_at_the_end = glbl_params.kinematics_at_the_end
        self.chunk_size = glbl_params.chunk_size
        self.ref_band = glbl_params.maglim[0]
        self.position = position
        self.max_mass = max_mass
        self.logger = logger

        # SPISEA specific setings
        self.spisea_dir = spisea_dir
        os.mkdirs(self.spisea_dir, exist_ok=True)
        self.spisea_evo_model = getattr(spisea.evolution, evolution_module.evo_model_name)()
        self.spisea_atm_func = getattr(spisea.atmospheres, evolution_module.atm_func_name)()
        self.spisea_wd_atm_func = getattr(spisea.atmospheres, evolution_module.wd_atm_func_name)()
        self.spisea_ifmr = None
        if evolution_module.ifmr_name is not None:
            self.spisea_ifmr = getattr(spisea.ifmr, evolution_module.ifmr_name)()
        self.spisea_imf = imf_module.spisea_imf
        self.multiplicity = imf_module.add_companions
        self.spisea_imf = getattr(spisea.imf.imf, imf_module.spisea_imf_name)()
        self.spisea_filters = glbl_params.chosen_bands
        
        # SPISEA limitations
        self.spisea_max_log_age = 10.14
        self.spisea_min_log_age = 5.0101
        self.spisea_age_bins =  100

    def generate_stars(self, radii, missing_stars, mass_limit,
        do_kinematics, props, avg_mass_per_star=None):
        position = np.vstack([
            np.column_stack(self.position.draw_random_point_in_slice(r_inner, r_outer, n_stars))
            for r_inner, r_outer, n_stars in zip(radii, radii[1:], missing_stars)
            ])

        min_mass = min(mass_limit)
        r_inner = np.repeat(radii[:-1], missing_stars)

        if self.kinematics_at_the_end:
            proper_motions = np.full((len(position), 3), np.nan)
            velocities = np.full((len(position), 3), np.nan)
            vr_lsr = np.repeat(np.nan, len(position))
        else:
            u, v, w, vr_hc, mu_l, mu_b, vr_lsr = do_kinematics(
                position[:, 3], position[:, 4], position[:, 5],
                position[:, 0], position[:, 1], position[:, 2]
                )
            proper_motions = np.column_stack([vr_hc, mu_l, mu_b])
            velocities = np.column_stack([u, v, w, ])

        # generate star at the positions
        return position, r_inner, proper_motions, velocities, vr_lsr, \
            self.generate_star_at_location(
            position[:, 0:3], props, min_mass, self.max_mass, avg_mass_per_star=avg_mass_per_star)

    def generate_star_at_location(self, position, props, min_mass=None, max_mass=None, avg_mass_per_star=None):
        """
        generates stars at the given positions
        """
        n_stars = len(position)

        # First - check whether the age distribution is uniform
        single_age = (self.age_module.age_func_name=='single_value')
        single_feh = (self.met_module.metallicity_func_name=='single_value')
        
        # NOTE: MAY NEED TO DEAL WITH FE/H vs M/H
        clusters=[]
        if single_age and single_feh:
            # First - determine whether to recomp
            isochrone = spisea.synthetic.IsochronePhot(logAge=np.log10(self.age_module.age_value*1e9), AKs=0,
                                distance=10, metallicity=self.met_module.metallicity_value,
                                evo_model=self.evo_model, atm_func=self.spisea_atm_func,
                                wd_atm_func=self.spisea_wd_atm_func, iso_dir=self.spisea_dir,
                                min_mass=min_mass, max_mass=max_mass,
                                filters=self.spisea_filters)
            cluster = spisea.synthetic.ResolvedCluster(isochrone, self.spisea_imf, n_stars*avg_mass_per_star, ifmr=self.spisea_ifmr)
            clusters.append(cluster)
        elif single_age:
            age = self.age_module.age_value
            fehs = self.met_module.draw_random_metallicity(
                    N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=age)
            pass
        elif single_metallicity:
            ages = self.age_module.draw_random_age(n_stars)
            feh = self.met_module.metallicity_value
            ages_hist, ages_bin_edges = np.histogram(ages, bins=self.spisea_age_bins,
                                        range=(10**self.spisea_min_log_age/1e9, 10**self.spisea_max_log_age/1e9))
            ages_hist_run = ages_hist[ages_hist>0]
            ages_bin_vals = np.diff(ages_bin_edges)[ages_hist>0]
            for i,age_bin in enumerate(ages_bin_vals):
                isochrone = spisea.synthetic.IsochronePhot(logAge=np.log10(self.age_module.age_value*1e9), AKs=0,
                                distance=10, metallicity=self.met_module.metallicity_value,
                                evo_model=self.evo_model, atm_func=self.spisea_atm_func,
                                wd_atm_func=self.spisea_wd_atm_func, iso_dir=self.spisea_dir,
                                min_mass=min_mass, max_mass=max_mass,
                                filters=self.spisea_filters)
                cluster = spisea.synthetic.ResolvedCluster(isochrone, self.spisea_imf, ages_hist_run[i]*avg_mass_per_star, ifmr=self.spisea_ifmr)
                clusters.append(cluster)
        else:
            ages = self.age_module.draw_random_age(n_stars)
            fehs = self.met_module.draw_random_metallicity(
                N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=ages)
            pass
            

        ref_mag, s_props, final_phase_flag, inside_grid, not_evolved = self.get_evolved_props(
            m_initial, met, age, props)

        return m_initial, age, met, ref_mag, s_props, final_phase_flag, inside_grid, not_evolved

