
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
from spisea.imf import imf as spisea_imf
from spisea.imf import multiplicity as spisea_multiplicity
from spisea import synthetic as spisea_synthetic
from spisea import ifmr as spisea_ifmr
import os, sys
import pdb
import pandas as pd
from astropy.table import vstack

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

class BlockSpiseaPrints:
    def __init__(self, block_prints):
        self.block_prints=block_prints

    def __enter__(self):
        if self.block_prints:
            self._original_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.block_prints:
            sys.stdout.close()
            sys.stdout = self._original_stdout

class SpiseaGenerator(StarGenerator):
    def __init__(self, density_module, imf_module, age_module, met_module, evolution_module,
            glbl_params, position, max_mass, ifmr_module, mult_module, bands, logger):
        # General synthpop things
        spisea_dir=const.ISOCHRONES_DIR+'/spisea/'
        self.generator_name = 'SpiseaGenerator'
        self.density_module = density_module
        self.imf_module = imf_module
        self.ifmr_module = ifmr_module
        self.mult_module = mult_module
        self.age_module = age_module
        self.met_module = met_module
        self.evolution_module = evolution_module
        self.kinematics_at_the_end = glbl_params.kinematics_at_the_end
        self.chunk_size = glbl_params.chunk_size
        self.ref_band = glbl_params.maglim[0]
        self.bands = bands
        self.obsmag = glbl_params.obsmag
        self.position = position
        self.max_mass = max_mass
        self.logger = logger
        self.system_mags = True

        # SPISEA specific setings
        self.spisea_dir = spisea_dir+evolution_module.spisea_evolution.model_version_name+'/'
        os.makedirs(self.spisea_dir, exist_ok=True)

        if self.evolution_module.name != 'SpiseaCluster':
            raise ValueError("To use SpiseaGenerator, the evolution class must be SpiseaCluster.")        

        if (self.mult_module is not None) and (self.mult_module.name!='SpiseaMultiplicity'):
            raise ValueError("Only SpiseaMultiplicity Multiplicity objects can be used by SpiseaGenerator")
        if self.imf_module.imf_name=='Kroupa':
            if self.mult_module is not None:
                self.imf_module.spisea_imf = spisea_imf.Kroupa_2001(multiplicity=self.mult_module.spisea_multiplicity)
            else:
                self.imf_module.spisea_imf = spisea_imf.Kroupa_2001()
        elif self.imf_module.imf_name=='PiecewisePowerlaw':
            if self.mult_module is not None:
                self.imf_module.spisea_imf = spisea_imf.IMF_broken_powerlaw([imf_module.min_mass, *imf_module.splitpoints, imf_module.max_mass],
                                                    -np.array(imf_module.alphas), multiplicity=self.mult_module.spisea_multiplicity)
            else:
                self.imf_module.spisea_imf = spisea_imf.IMF_broken_powerlaw([imf_module.min_mass, *imf_module.splitpoints, imf_module.max_mass],
                                                    -np.array(imf_module.alphas))
        elif self.imf_module.imf_name=='SpiseaImf':
            if (self.imf_module.spisea_multiplicity is None) and (self.mult_module is None):
                pass
            # Re-initialize IMF module with proper multiplicity if needed
            elif (self.imf_module.spisea_multiplicity is None) and (self.mult_module.multiplicity_name=='SpiseaMultiplicity'):
                self.imf_module.spisea_imf = getattr(spisea_imf, self.imf_module.spisea_imf_name)(
                            massLimits=np.array([self.imf_module.min_mass, self.imf_module.max_mass]), 
                            multiplicity=self.mult_module.spisea_multiplicity, **self.imf_module.spisea_imf_kwargs)
        else:
            raise ValueError("Invalid IMF for SPISEA Generator; must use Kroupa, PiecewisePowerlaw, or SpiseaImf")

        self.mh_list = np.log10(np.array(self.evolution_module.spisea_evolution.z_list) / self.evolution_module.spisea_evolution.z_solar)

    @staticmethod
    def spisea_props_to_synthpop(spisea_tab):
        renames = {'mass_current': 'Mass', 'mass':'iMass', 'logg':'log_g',
                          'N_companions': 'n_companions', 'e':"eccentricity"}
        sp_dict = {}
        for col in spisea_tab.columns:
            if col in renames:
                sp_dict[renames[col]] = np.array(spisea_tab[col])
            elif col not in ['Teff', 'L', 'isMultiple', 'systemMass']:
                sp_dict[col] = np.array(spisea_tab[col])
        lums = np.array(spisea_tab['L'])
        teffs = np.array(spisea_tab['Teff'])
        sp_dict['log_L'] = np.log10(lums/const.Lsun_w)
        sp_dict['log_Teff'] = np.log10(teffs)
        sp_dict['log_R'] = np.log10((np.sqrt(lums/(4*np.pi*const.sigma_sb*teffs**4)))/const.Rsun_m)
        sp_dict['star_mass'] = sp_dict['Mass']

        return sp_dict

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
                    sys_idxs = star_systems['system_idx'].to_numpy()
                    dist_modulus_series = pd.Series(dist_modulus, index=sys_idxs)
                    companions.loc[:,band] += dist_modulus_series[
                            companions['system_idx'].to_numpy()].to_numpy()

        return star_systems, companions

    def generate_star_at_location(self, position, props, min_mass=None, max_mass=None, avg_mass_per_star=None):
        """
        generates stars at the given positions
        """
        #pdb.set_trace()
        if avg_mass_per_star is None:
            avg_mass_per_star = 1 #self.synthpop_imf_module.average_mass(min_mass=min_mass, max_mass=self.max_mass)

        n_stars = len(position)

        # First - check whether the age distribution is uniform
        single_age = (self.age_module.age_func_name=='single_value')
        single_feh = (self.met_module.metallicity_func_name=='single_value')
        
        # NOTE: We check Fe/H then convert to M/H for SPISEA
        if single_age and single_feh:
            age_all = self.evolution_module.log_age_list[np.argmin(np.abs(self.evolution_module.log_age_list-np.log10(self.age_module.age_value*1e9)))]
            mh_all = self.mh_list[np.argmin(np.abs(self.evolution_module.feh_list-self.met_module.metallicity_value))]
            comb_bin_idxs = np.arange(n_stars)
            bins2d = [[age_all, feh_all, n_stars]]
        elif single_age:
            # Sample metallicities in [Fe/H], then bin by nearest grid point
            age_all = self.evolution_module.log_age_list[np.argmin(np.abs(self.evolution_module.log_age_list-np.log10(self.age_module.age_value*1e9)))]
            fehs = self.met_module.draw_random_metallicity(
                    N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=10**age_all/1e9)
            feh_bins, comb_bin_idxs, feh_bin_cts = np.unique(np.argmin(np.abs(self.evolution_module.feh_list - fehs[:, None]), axis=1), 
                                                              return_inverse=True, return_counts=True)
            bins2d = np.transpose([np.ones(len(feh_bins))*age_all, self.mh_list[feh_bins], feh_bin_cts])
        elif single_feh:
            # Sample ages in log10(yr), then bin by nearest grid point
            ages = np.log10(self.age_module.draw_random_age(n_stars)*1e9)
            age_bins, comb_bin_idxs, age_bin_cts = np.unique(np.argmin(np.abs(self.evolution_module.log_age_list - ages[:, None]), axis=1), 
                                                              return_inverse=True, return_counts=True)
            mh_all = self.mh_list[np.argmin(np.abs(self.evolution_module.feh_list-self.met_module.metallicity_value))]
            bins2d = np.transpose([self.evolution_module.log_age_list[age_bins], np.ones(len(age_bins))*mh_all, age_bin_cts])
        else:
            ages = np.log10(self.age_module.draw_random_age(n_stars)*1e9)
            fehs = self.met_module.draw_random_metallicity(
                N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=10**np.mean(ages)/1e9)
            comb_vals = np.transpose([np.argmin(np.abs(self.evolution_module.log_age_list - ages[:, None]), axis=1), np.argmin(np.abs(self.evolution_module.feh_list-fehs[:,None]), axis=1)])
            comb_bins, comb_bin_idxs, comb_bin_cts = np.unique(comb_vals, axis=0, return_inverse=True, return_counts=True)
            bin_ages = self.evolution_module.log_age_list[comb_bins[:,0]]
            bin_mhs = self.mh_list[comb_bins[:,1]]
            bins2d = np.transpose([bin_ages, bin_mhs, comb_bin_cts])

        # Set up data arrays
        star_systems_list = []
        companions_list = []
        # TODO: need to deal with matching stars back in to proper positions in case metallicity is position-dependent :/
        # We have the bin array inverses now, so that should help, right ??
        star_systems_data = {param: np.zeros(n_stars) for param in props}
        star_systems_data['iMass'] = np.zeros(n_stars)
        star_systems_data['Mass'] = np.zeros(n_stars)
        star_systems_data['age'] = np.zeros(n_stars)
        star_systems_data['Fe/H_initial'] = np.zeros(n_stars)
        star_systems_data['system_idx'] = np.zeros(n_stars, dtype=int)
        star_systems_data['n_companions'] = np.zeros(n_stars)
        companions_data = {param: [] for param in props}
        companions_data['system_idx'] = []
        companions_data['iMass'] = []
        companions_data['Mass'] = []
        companions_data['eccentricity'] = []
        companions_data['log_a'] = []

        max_system_idx = -1
        for i_bin, bin2d in enumerate(bins2d):
            # Figure out where stars in this bin fit in the data set (so their positions and age/met can be correlated)
            star_idxs_in_bin = np.where(comb_bin_idxs==i_bin)[0]
            # Use a minimum mass per cluster so we don't get an error
            star_systems_list_bin = []
            companions_list_bin = []
            n_bin = int(bin2d[-1])
            self.logger.debug("Starting SPISEA cluster generation for bin log_age="+str(round(bin2d[0],2))+
                                " [M/H]="+str(bin2d[1])+" for "+str(n_bin)+" stars")
            cluster_stars_needed = n_bin
            generate_mass = np.minimum(np.maximum(cluster_stars_needed*avg_mass_per_star*1.1, 100.0), 
                                    self.chunk_size*avg_mass_per_star)
            # Loop until we have enough stars
            while cluster_stars_needed > 0:
                with BlockSpiseaPrints(block_prints=self.evolution_module.block_spisea_prints):
                    isochrone = spisea_synthetic.IsochronePhot(logAge=bin2d[0], AKs=0,
                                        distance=10, metallicity=bin2d[1],
                                        evo_model=self.evolution_module.spisea_evolution, atm_func=self.evolution_module.spisea_atm_func,
                                        wd_atm_func=self.evolution_module.spisea_wd_atm_func, iso_dir=self.spisea_dir,
                                        min_mass=np.min(min_mass), max_mass=max_mass,
                                        filters=self.evolution_module.bands_obs_str)
                    cluster=spisea_synthetic.ResolvedCluster(isochrone, self.imf_module.spisea_imf, generate_mass, 
                                                    ifmr=self.ifmr_module.spisea_ifmr, keep_low_mass_stars=True)
                star_systems_i = cluster.star_systems
                if "companions" in cluster.__dir__():
                    companions_i = cluster.companions
                    companions_i['system_idx'] += (max_system_idx + 1)
                    companions_list_bin.append(companions_i)
                else:
                    star_systems_i['n_companions'] = 0
                star_systems_i['system_idx'] = np.arange(len(star_systems_i)) + max_system_idx + 1
                max_system_idx = star_systems_i['system_idx'].max()
                keep_idx = ((star_systems_i['mass']>np.min(min_mass)) & (star_systems_i['mass']<max_mass))
                star_systems_i = star_systems_i[keep_idx]
                star_systems_list_bin.append(star_systems_i)
                cluster_stars_needed -= len(star_systems_i)

            star_systems_bin = vstack(star_systems_list_bin)
            if len(companions_list_bin)>0:
                companions_bin = vstack(companions_list_bin)
            else:
                companions_bin = None
            # Drop any excess stars
            if cluster_stars_needed<0:
                star_systems_bin = star_systems_bin[:cluster_stars_needed]
            # Get the data into the expected form and dropped in place in the star list
            star_systems_bin = self.spisea_props_to_synthpop(star_systems_bin)
            for param in list(props)+['iMass','Mass','system_idx', 'n_companions']:
                star_systems_data[param][star_idxs_in_bin] = star_systems_bin[param]
            star_systems_data['age'][star_idxs_in_bin] = 10**bin2d[0] / 1e9
            star_systems_data['Fe/H_initial'][star_idxs_in_bin] = self.evolution_module.feh_list[np.argmin(np.abs(self.mh_list-bin2d[1]))]
            # Drop any companions whose systems got dropped
            if self.imf_module.spisea_imf.make_multiples and len(companions_bin)>0:
                companions_bin = companions_bin[np.isin(companions_bin['system_idx'], star_systems_data['system_idx'])]
                companions_bin = self.spisea_props_to_synthpop(companions_bin)
            if self.imf_module.spisea_imf.make_multiples and len(companions_bin)>0:
                for param in list(props)+['iMass','Mass','system_idx', 'eccentricity', 'log_a']:
                    companions_data[param] += list(companions_bin[param])

        star_systems = pd.DataFrame(star_systems_data)
        star_systems.drop(columns='star_mass', inplace=True)
        companions = pd.DataFrame(companions_data) if (self.imf_module.spisea_imf.make_multiples) else None

        star_systems.loc[:,'system_Mass'] = star_systems['Mass']
        if self.imf_module.spisea_imf.make_multiples:
            #companions = self.spisea_props_to_synthpop(companions)
            companions.drop(columns='star_mass', inplace=True)
            if len(companions)>0:
                comp_mass_sums = companions.groupby("system_idx")['Mass'].sum()
                primary_idxs = star_systems.index[star_systems['n_companions']>0]
                star_systems.loc[primary_idxs,'system_Mass'] += comp_mass_sums[
                                        star_systems['system_idx'][primary_idxs]]

        return star_systems, companions

