
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
from spisea import evolution as spisea_evolution
from spisea import atmospheres as spisea_atmospheres
from spisea.imf import imf as spisea_imfs
from spisea import synthetic as spisea_synthetic
from spisea import ifmr as spisea_ifmr
import os
import pdb
import pandas as pd

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
    def __init__(self, density_module, imf_module, age_module, met_module, evolution_module,
            glbl_params, position, max_mass, ifmr_module, mult_module, bands, logger):
        # General synthpop things
        spisea_dir=const.ISOCHRONES_DIR+'/spisea/'
        self.generator_name = 'SpiseaGenerator'
        self.density_module = density_module
        self.age_module = age_module
        self.met_module = met_module
        self.ifmr_module = ifmr_module
        self.mult_module = mult_module
        self.kinematics_at_the_end = glbl_params.kinematics_at_the_end
        self.chunk_size = glbl_params.chunk_size
        self.ref_band = glbl_params.maglim[0]
        self.position = position
        self.max_mass = max_mass
        self.logger = logger
        self.bands = bands
        self.system_mags = True
        self.synthpop_imf_module = imf_module

        # SPISEA specific setings
        if evolution_module.evo_model_name[:6]=='MISTv1': #.2') or (evolution_module.evo_model_name=='MISTv1.0'):
            n_version = float(evolution_module.evo_model_name[5:8])
            synthpop_extension = ("synthpop" in evolution_module.evo_model_name.lower())
            self.spisea_evo_model = spisea_evolution.MISTv1(version=n_version,
                                                            synthpop_extension=synthpop_extension)
            self.feh_list = np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.75,-1.5,-1.25,
                                      -1.0,-0.75,-0.5,-0.25,0,0.25,0.5])
            self.log_age_list = np.linspace(5.0,10.3,107)
            self.log_age_list[0] = 5.01
        elif evolution_module.evo_model_name=='MergedBaraffePisaEkstromParsec':
            self.spisea_evo_model = spisea_evolution.MergedBaraffePisaEkstromParsec()
            self.feh_list = np.array([0.0])
            self.log_age_list = np.linspace(6.0,10.1,83)
            self.log_age_list[-1] = 10.9
            Warning(f"Evolution module {evolution_module.evo_model_name} only includes solar metallicy. All stars will be assigned solar metallicity.")
        else:
            raise ValueError("Invalid SPISEA evolution_model. Only MISTv1.0, MISTv1.2, and MergedBaraffePisaEkstromParsec are available at this time.")
        self.spisea_dir = spisea_dir+evolution_module.evo_model_name+'/'
        os.makedirs(self.spisea_dir, exist_ok=True)

        #self.spisea_evo_model = getattr(spisea_evolution, evolution_module.evo_model_name)()
        self.spisea_atm_func = getattr(spisea_atmospheres, evolution_module.atm_func_name)
        self.spisea_wd_atm_func = getattr(spisea_atmospheres, evolution_module.wd_atm_func_name)
        #self.spisea_ifmr = None
        if self.ifmr_module.name == 'PopsycleIfmrs':
            self.spisea_ifmr = getattr(spisea_ifmr, ifmr_module.spisea_ifmr_name)()
        else:
            raise ValueError("Invalid IFMR module for SpiseaGenerator. Only PopsycleIfmrs is available at this time.")
        self.spisea_multiplicity = None
        if evolution_module.multiplicity_name is not None:
            self.spisea_multiplicity = getattr(spisea.multiplicity, evolution_module.multiplicity_name)()
            raise NotImplementedError("Stellar multiplicity via SPISEA is not yet implemented.")
        if imf_module.imf_name=='Kroupa':
            self.spisea_imf = spisea_imfs.Kroupa_2001(multiplicity=self.spisea_multiplicity)
        elif imf_module.imf_name=='Piecewise Powerlaw':
            self.spisea_imf = spisea_imfs.IMF_broken_powerlaw([imf_module.min_mass, *imf_module.splitpoints, imf_module.max_mass],
                                                    -np.array(imf_module.alphas), multiplicity=self.spisea_multiplicity)
        else:
            raise ValueError("Invalid IMF for SPISEA Generator; must use Kroupa or PiecewisePowerlaw.")
        #self.spisea_imf = imf_module.spisea_imf
        #self.multiplicity = imf_module.add_companions
        #self.spisea_imf = getattr(spisea_imf.imf, imf_module.spisea_imf_name)()
        self.spisea_filters = [filt.replace('-',',') for filt in evolution_module.bands]
        
        # SPISEA limitations
        self.spisea_max_log_age = 10.14
        self.spisea_min_log_age = 5.0101
        self.spisea_max_feh = 0.5
        self.spisea_min_feh = -4.0
        #self.spisea_age_bins =  100
        self.mh_list = np.log10(np.array(self.spisea_evo_model.z_list) / self.spisea_evo_model.z_solar)

    @staticmethod
    def spisea_props_to_synthpop(spisea_df):
        #print(spisea_df.keys())
        spisea_df.rename({'mass_current': 'Mass', 'mass':'iMass', 'logg':'log_g',
                          'N_companions': 'n_companions'})
        pd.eval('log_L = np.log10(lums/Lsun_w)', local_dict={'Lsun_w': const.Lsun_w,
                'lums':spisea_df['L']}, target=spisea_df, inplace=True)
        pd.eval('log_Teff = np.log10(Teffs)', local_dict={'Teffs'=spisea_df['Teff']},
                target=spisea_df, inplace=True)
        pd.eval('log_R = np.log10((np.sqrt(lums/(4*np.pi*sigma_sb*Teffs**4)))/Rsun_m)',
                local_dict={'Lsun_w': const.Lsun_w, 'sigma_sb':const.sigma_sb, 'Rsun_m':const.Rsun_m,
                'lums':spisea_df['L'], 'Teffs'=spisea_df['Teff']},
                target=spisea_df, inplace=True)
        spisea_df.drop(columns=['Teff', 'L', 'isMultiple', 'systemMass'], inplace=True, errors='ignore')

        return spisea_df

    def generate_stars(self, radii, missing_stars, mass_limit,
        do_kinematics, props):
        position = np.vstack([
            np.column_stack(self.position.draw_random_point_in_slice(r_inner, r_outer, n_stars, population_density_func=self.density_module.density))
            for r_inner, r_outer, n_stars in zip(radii, radii[1:], missing_stars)
            ])

        min_mass = np.repeat(mass_limit, missing_stars)
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
        star_systems, companions = self.generate_star_at_location(
            position[:, 0:3], props, min_mass, self.max_mass)
        star_systems.loc[:,'x'] = position[:,0]
        star_systems.loc[:,'y'] = position[:,1]
        star_systems.loc[:,'z'] = position[:,2]
        star_systems.loc[:,'Dist'] = position[:,3]
        star_systems.loc[:,'l'] = position[:,4]
        star_systems.loc[:,'b'] = position[:,5]
        star_systems.loc[:,'r_inner'] = r_inner
        star_systems.loc[:,'vr_bc'] = proper_motions[:,0]
        star_systems.loc[:,'mul'] = proper_motions[:,1]
        star_systems.loc[:,'mub'] = proper_motions[:,2]
        star_systems.loc[:,'U'] = velocities[:,0]
        star_systems.loc[:,'V'] = velocities[:,1]
        star_systems.loc[:,'W'] = velocities[:,2]
        star_systems.loc[:,'VR_LSR'] = vr_lsr
        
        if self.obsmag:
            dist_modulus = 5*np.log10(position[:,3] * 100)
            for band in self.bands:
                star_systems.loc[:,band] += dist_modulus
                if companions is not None:
                    companions.loc[:,band] += dist_modulus[
                            companions['system_idx'].to_numpy()]

        return star_systems, companions

    def generate_star_at_location(self, position, props, min_mass=None, max_mass=None, avg_mass_per_star=None):
        """
        generates stars at the given positions
        """
        #pdb.set_trace()
        if avg_mass_per_star is None:
            avg_mass_per_star = self.synthpop_imf_module.average_mass(min_mass=min_mass, max_mass=self.max_mass)

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
                    N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=10**age_all/1e9)
            feh_bins = np.unique(np.argmin(np.abs(self.feh_list - fehs[:, None]), axis=1), return_counts=True)
            bins2d = np.transpose([np.ones(len(feh_bins[0]))*age_all, self.mh_list[feh_bins[0]], feh_bins[1]])
        elif single_feh:
            # Sample ages in log10(yr), then bin by nearest grid point
            ages = np.log10(self.age_module.draw_random_age(n_stars)*1e9)
            age_bins = np.unique(np.argmin(np.abs(self.log_age_list - ages[:, None]), axis=1), return_counts=True)
            mh_all = self.mh_list[np.argmin(np.abs(self.feh_list-self.met_module.metallicity_value))]
            bins2d = np.transpose([self.log_age_list[age_bins[0]], np.ones(len(age_bins[0]))*mh_all, age_bins[1]])
        else:
            ages = np.log10(self.age_module.draw_random_age(n_stars)*1e9)
            fehs = self.met_module.draw_random_metallicity(
                N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=10**np.mean(ages)/1e9)
            comb_vals = np.transpose([np.argmin(np.abs(self.log_age_list - ages[:, None]), axis=1), np.argmin(np.abs(self.feh_list-fehs[:,None]), axis=1)])
            comb_bins = np.unique(comb_vals, axis=0, return_counts=True)
            bin_ages = self.log_age_list[comb_bins[0][:,0]]
            bin_mhs = self.mh_list[comb_bins[0][:,1]]
            bins2d = np.transpose([bin_ages, bin_mhs, comb_bins[1]])

        # Set up data arrays
        star_systems_list = []
        companions_list = []
        # TODO: need to deal with matching stars back in to proper positions in case metallicity is position-dependent :/

        max_system_idx = -1
        for bin2d in bins2d:
            # Use a minimum mass per cluster so we don't get an error
            star_systems_list_bin = []
            companions_list_bin = []
            n_bin = int(bin2d[-1])
            print("Starting SPISEA cluster generation for bin log_age="+str(round(bin2d[0],2))+" [M/H]="+str(bin2d[1])+" for "+str(n_bin)+" stars")
            cluster_stars_needed = n_bin
            generate_mass = np.minimum(np.maximum(cluster_stars_needed*avg_mass_per_star*1.1, 10.0), 
                                    self.chunk_size*avg_mass_per_star)
            # Loop until we have enough stars
            while cluster_stars_needed > 0:
                isochrone = spisea_synthetic.IsochronePhot(logAge=bin2d[0], AKs=0,
                                    distance=10, metallicity=bin2d[1],
                                    evo_model=self.spisea_evo_model, atm_func=self.spisea_atm_func,
                                    wd_atm_func=self.spisea_wd_atm_func, iso_dir=self.spisea_dir,
                                    min_mass=min_mass, max_mass=max_mass,
                                    filters=self.spisea_filters)
                cluster=spisea_synthetic.ResolvedCluster(isochrone, self.spisea_imf, generate_mass, ifmr=self.spisea_ifmr,
                                                            keep_low_mass_stars=True)
                star_systems_i = cluster.star_systems.to_pandas()
                if cluster.companions is not None:
                    companions_i = cluster.companions.to_pandas()
                keep_idx = ((star_systems_i['mass']>min_mass) & (star_systems_i['mass']<max_mass))
                star_systems_i = star_systems_i[keep_idx]
                if cluster.companions is not None:
                    companions_i.loc[:,'system_idx'] += (max_system_idx + 1)
                    companions_list_bin.append(companions_i)
                star_systems_i.loc[:, 'system_idx'] = np.arange(len(star_systems_i)) + max_system_idx + 1
                star_systems_list_bin.append(star_systems_i)
                max_system_idx = star_systems_i['system_idx'].max()
                cluster_stars_needed -= len(star_systems_i)

            star_systems_bin = pd.concat(clusters, ignore_index=True)
            if len(companions_list_bin)>0:
                companions_bin = pd.concat(companions_list_bin)
            else:
                companions_bin = None
            # Drop any excess stars
            if cluster_stars_needed<0:
                drop_idx = np.random.choice(star_systems_bin.index.to_numpy(),
                                size=-cluster_stars_needed, replace=False)
                star_systems_bin.drop(index=drop_idx, inplace=True)
            # Drop any companions whose systems got dropped
            if self.mult_module is not None and len(companions>0):
                companions_bin = companions_bin[np.isin(companions['system_idx'], star_systems['system_idx'])]
                
            star_systems.loc[:,'age'] = 10**bin2d[0] / 1e9
            star_systems.loc[:,'Fe/H_initial'] =self.feh_list[np.argmin(np.abs(self.mh_list-bin2d[1]))]

            # Bin complete
            star_systems_list.append(star_systems_bin)
            companions_list.append(companions_bin)
            
        star_systems = pd.concat(star_systems_list)
        companions = pd.concat(companions_list) if (self.mult_module is not None) else None
        
        star_systems = self.spisea_props_to_synthpop(star_systems)
        star_systems.loc[:,'system_Mass'] = star_systems['Mass']
        if self.mult_module is not None and len(companions>0):
            companions = self.spisea_props_to_synthpop(companions)
            comp_mass_sums = companions.groupby("system_idx")['Mass'].sum()
            primary_idxs = star_systems.index[star_systems['n_companions']>0]
            star_systems.loc[primary_idxs,'system_Mass'] += comp_mass_sums[
                                    star_systems['system_idx'][primary_idxs]]

        return star_systems, companions

