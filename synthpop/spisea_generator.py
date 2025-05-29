
"""
This file contains the in-progress Spisea alternative to StarGenerator,
It bins stars by age and metallicity and generates SPISEA clusters,
    then assigns them locations based on the population_density
"""

__all__ = ["SpiseaGenerator"]
__author__ = "M.J. Huston"
__credits__ = ["M.J. Huston"]
__date__ = "2025-05-28"

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
        self.generator_name = 'SpiseaGenerator'
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
        if (evolution_module.evo_model_name=='MISTv1.2') or (evolution_module.evo_model_name=='MISTv1.0'):
            self.evolution_model = evolution.MISTv1(version=float(evolution_model_name[-3:]))
            self.feh_list = np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.25,0.5])
            self.log_age_list = np.linspace(5.0,10.3,107)
            self.log_age_list[0] = 5.01
        else:
            raise ValueError("Invalid SPISEA evolution_model. Only MISTv1.0 and MISTv1.2 are available at this time.")

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
        self.spisea_max_feh = 0.5
        self.spisea_min_feh = -4.0
        #self.spisea_age_bins =  100
        self.mh_list = np.log10(self.spisea_evo_model.z_list / self.spisea_evo_model.z_solar)

    @staticmethod
    def spisea_props_to_synthpop(synth_props, spisea_df):
        conv_props = {}
        for prop in synth_props:
            if prop=='phase':
                conv_props[prop] = spisea_df['phase'].to_numpy()
            elif prop=='log_L':
                conv_props[prop] = np.log10(spisea_df['L'].to_numpy()/const.Lsun_w)
            elif prop=='log_Teff':
                conv_props[prop] = np.log10(spisea_df['Teff'].to_numpy())
            elif prop=='log_g':
                conv_props[prop] = spisea_df['logg'].to_numpy()
            elif prop'log_R':
                conv_props[prop] = np.log10(spisea_df['R'].to_numpy()/const.Rsun_m)
            else:
                try:
                    spisea_col = 'm_'+spisea.synthetic.get_filter_col_name(prop.replace('-',','))
                    conv_props[prop] = spisea_df[spisea_col].to_numpy()
                except:
                    raise ValueError("Invalid column for SPISEA generator: "+str(prop))

        return conv_props

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
        
        # NOTE: We check Fe/H then convert to M/H for SPISEA
        if single_age and single_feh:
            age_all = self.log_age_list[np.argmin(np.abs(self.log_age_list-np.log10(self.age_module.age_value*1e9)))]
            mh_all = self.mh_list[np.argmin(np.abs(self.feh_list-self.met_module.metallicity_value))]
            bins2d = [[age_all, feh_all, n_stars]]
        elif single_age:
            # Sample metallicities in [Fe/H], then bin by nearest grid point
            age_all = self.log_age_list[np.argmin(np.abs(self.log_age_list-np.log10(self.age_module.age_value*1e9)))]
            fehs = self.met_module.draw_random_metallicity(
                    N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=age)
            feh_bins = np.unique(np.argmin(np.abs(self.feh_list - fehs[:, None]), axis=1), return_counts=True)
            bins2d = np.transpose([np.ones(len(feh_bins[0]))*age_all, self.mh_list[feh_bins[0]], feh_bins[1]])
        elif single_metallicity:
            # Sample ages in log10(yr), then bin by nearest grid point
            ages = np.log10(self.age_module.draw_random_age(n_stars)*1e9)
            age_bins = np.unique(np.argmin(np.abs(self.log_age_list - ages[:, None]), axis=1), return_counts=True)
            mh_all = self.mh_list[np.argmin(np.abs(self.feh_list-self.met_module.metallicity_value))]
            bins2d = np.transpose([self.log_age_list[age_bins[0]], np.ones(len(age_bins[0]))*mh_all, age_bins[1]])
        else:
            ages = np.log10(self.age_module.draw_random_age(n_stars)*1e9)
            fehs = self.met_module.draw_random_metallicity(
                N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=ages)
            comb_vals = np.transpose([np.argmin(np.abs(self.log_age_list - ages[:, None]), axis=1), np.argmin(np.abs(self.feh_list-fehs[:,None]), axis=1)])
            comb_bins = np.unique(comb_vals, return_counts=True)
            bin_ages = self.log_age_list[comb_bins[0][:,0]]
            bin_mhs = self.mh_list[comb_bins[0][:,1]]
            bins2d = np.transpose([bin_ages, bin_mhs, comb_bins[1]])

        # Set up data arrays
        m_initial = np.zeros(n_stars)
        age = np.zeros(n_stars)
        met = np.zeros(n_stars)
        ref_mag = np.zeros(n_stars)
        s_props = {p: np.zeros(n_stars) for p in props}
        final_phase_flag = np.zeros(n_stars, bool)
        inside_grid = np.ones(n_stars, bool)
        not_evolved = np.zeros(n_stars, bool)

        stars_done = 0
        for bin2d in bins2d:
            # Use a minimum mass per cluster so we don't get an error
            clusters = []
            cluster_stars_needed = bin2d[2]
            generate_mass = np.min(np.max(cluster_stars_needed*avg_mass_per_star*1.1, 10.0), 
                                    self.chunk_size*avg_mass_per_star)
            # Loop until we have enough stars
            while cluster_stars_needed > 0:
                isochrone = spisea.synthetic.IsochronePhot(logAge=bin2d[0], AKs=0,
                                    distance=10, metallicity=bin2d[1],
                                    evo_model=self.evo_model, atm_func=self.spisea_atm_func,
                                    wd_atm_func=self.spisea_wd_atm_func, iso_dir=self.spisea_dir,
                                    min_mass=min_mass, max_mass=max_mass,
                                    filters=self.spisea_filters)
                clusters.append(spisea.synthetic.ResolvedCluster(isochrone, self.spisea_imf, mass_needed, ifmr=self.spisea_ifmr))
                assert (not hasattr(clusters[-1], 'companions')), "Error: Companions not yet implemented."
                cluster_stars_needed -= len(cluster)
            cluster_comb = pd.concat(clusters, ignore_index=True)
            # Drop any excess stars
            if cluster_stars_needed<0:
                drop_idx = np.random.choice(cluster_comb.index.to_numpy(), size=-cluster_stars_needed)
                cluster_comb.drop(index=drop_idx, inplace=True)
            # Get data from SPISEA cluster into SynthPop's formats
            m_initial[stars_done:stars_done+bin2d[2]] = cluster_comb['mass'].to_numpy()
            age[stars_done:stars_done+bin2d[2]] = 10**bin2d[0] / 1e9
            met[stars_done:stars_done+bin2d[2]] = self.feh_list[np.argmin(np.abs(self.mh_list-bin2d[1]))]
            conv_props = self.spisea_props_to_synthpop(props)
            for prop in props:
                s_props[prop][stars_done:stars_done+bin2d[2]] = conv_props[prop]
            not_evolved[stars_done:stars_done+bin2d[2]] = np.isnan(cluster_comb.L.to_numpy())
            # Bin complete
            stars_done += bin2d[2]

        ref_mag = s_props[self.ref_band]
        return m_initial, age, met, ref_mag, s_props, final_phase_flag, inside_grid, not_evolved

