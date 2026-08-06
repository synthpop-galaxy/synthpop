"""
This file contains the SPISEA-based alternative to StarGenerator. 
It bins stars by age and metallicity and generates SPISEA clusters,
    then assigns them locations based on the population_density.
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
from multiprocessing import Pool

# Local Imports
# used to allow running as main and importing to another script
try:
    from . import constants as const
    from . import synthpop_utils as sp_utils
    from .star_generator import StarGenerator
    from .synthpop_utils.synthpop_logging import logger
    from .synthpop_utils import coordinates_transformation as coord_trans
    from .synthpop_utils import Parameters
    from .modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from .modules.age import Age
    from .modules.initial_mass_function import InitialMassFunction
    from .modules.kinematics import Kinematics
    from .modules.metallicity import Metallicity
    from .modules.population_density import PopulationDensity
except ImportError:
    import constants as const
    import synthpop_utils as sp_utils
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
            glbl_params, max_mass, ifmr_module, mult_module, bands, logger):
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
        self.chunk_size = glbl_params.chunk_size
        self.bands = bands
        self.obsmag = glbl_params.obsmag
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

        # TODO clean this up later
        if self.evolution_module.spisea_evolution.model_version_name=='COSMIC':
            self.mh_list = self.evolution_module.feh_list
        else:
            self.mh_list = np.log10(np.array(self.evolution_module.spisea_evolution.z_list) / self.evolution_module.spisea_evolution.z_solar)
        
        self.n_proc = evolution_module.n_proc

    def generate_star_at_location(self, position, props, 
                min_mass=None, max_mass=None, radii=None, avg_mass_per_star=None, skip_lowmass_stars=None):
        """
        Generates stars at the given positions
        """
        if avg_mass_per_star is None:
            avg_mass_per_star = 1 #self.synthpop_imf_module.average_mass(min_mass=min_mass, max_mass=self.max_mass)

        n_stars = len(position[0])

        # First - check whether the age distribution is uniform
        single_age = (self.age_module.age_func_name=='single_value')
        single_feh = (self.met_module.metallicity_func_name=='single_value')
        
        # NOTE: We check Fe/H then convert to M/H for SPISEA
        if single_age and single_feh:
            age_all = self.evolution_module.log_age_list[np.argmin(np.abs(self.evolution_module.log_age_list-np.log10(self.age_module.age_value*1e9)))]
            mh_all = self.mh_list[np.argmin(np.abs(self.evolution_module.feh_list-self.met_module.metallicity_value))]
            comb_bin_idxs = np.zeros(n_stars,dtype=int)
            bins2d = [[age_all, mh_all, n_stars]]
        elif single_age:
            # Sample metallicities in [Fe/H], then bin by nearest grid point
            age_all = self.evolution_module.log_age_list[np.argmin(np.abs(self.evolution_module.log_age_list-np.log10(self.age_module.age_value*1e9)))]
            fehs = self.met_module.draw_random_metallicity(
                    N=n_stars, x=position[0], y=position[1], z=position[2], age=10**age_all/1e9)
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
                N=n_stars, x=position[0], y=position[1], z=position[2], age=10**np.mean(ages)/1e9)
            comb_vals = np.transpose([np.argmin(np.abs(self.evolution_module.log_age_list - ages[:, None]), axis=1), np.argmin(np.abs(self.evolution_module.feh_list-fehs[:,None]), axis=1)])
            comb_bins, comb_bin_idxs, comb_bin_cts = np.unique(comb_vals, axis=0, return_inverse=True, return_counts=True)
            bin_ages = self.evolution_module.log_age_list[comb_bins[:,0]]
            bin_mhs = self.mh_list[comb_bins[:,1]]
            bins2d = np.transpose([bin_ages, bin_mhs, comb_bin_cts])
        
        # Serial case
        if self.n_proc==1:
            res = []
            for i_bin, bin2d in enumerate(bins2d):
                system_idxs = np.where(comb_bin_idxs==i_bin)[0]
                res.append(generate_spisea_cluster_stars(system_idxs, bin2d[0], bin2d[1],
                                       self.evolution_module.feh_list[np.argmin(np.abs(self.mh_list-bin2d[1]))],
                                       avg_mass_per_star,
                                       self.evolution_module.block_spisea_prints,
                                       self.evolution_module.spisea_evolution,
                                       self.evolution_module.spisea_atm_func,
                                       self.evolution_module.spisea_wd_atm_func,
                                       np.min(min_mass), max_mass,
                                       self.evolution_module.bands_obs_str,
                                       self.imf_module.spisea_imf,
                                       self.ifmr_module.spisea_ifmr,
                                       self.spisea_dir, props,
                                       self.evolution_module.bands,
                                       self.evolution_module.bbh_frac))                       
        
        # Parallel case:
        else:
            # Split up some bins for efficiency if relevant...
            if np.any(bins2d[:,2]>self.chunk_size):
                bins2d_new = []
                comb_bin_idxs_new = []
                for i_bin, bin2d in enumerate(bins2d):
                    comb_bin_idxs_split = np.split(comb_bin_idxs[i_bin],
                                            int(np.ciel(self.chunk_size/bin2d[2])))
                    bins2d_new += [[bin2d[0],bin2d[1],len(idxs)] for idxs in comb_bin_idxs_split]
                    comb_bin_idxs_new += comb_bin_idxs_split
                bins2d = bins2d_new
                comb_bin_idxs = comb_bin_idxs_new
                
            # Iterate through the bins and set up inputs for generate_spisea_cluster_stars
            cluster_inputs = []
            for i_bin, bin2d in enumerate(bins2d):
                system_idxs = np.where(comb_bin_idxs==i_bin)[0]
                cluster_inputs.append((system_idxs, bin2d[0], bin2d[1],
                                       self.evolution_module.feh_list[np.argmin(np.abs(self.mh_list-bin2d[1]))],
                                       avg_mass_per_star,
                                       self.evolution_module.block_spisea_prints,
                                       self.evolution_module.spisea_evolution,
                                       self.evolution_module.spisea_atm_func,
                                       self.evolution_module.spisea_wd_atm_func,
                                       np.min(min_mass), max_mass,
                                       self.evolution_module.bands_obs_str,
                                       self.imf_module.spisea_imf,
                                       self.ifmr_module.spisea_ifmr,
                                       self.spisea_dir, props,
                                       self.evolution_module.bands,
                                       self.evolution_module.bbh_frac))
                                       
            # Do the parallel generation
            with Pool(self.n_proc) as p:
                res = p.starmap_async(generate_spisea_cluster_stars, cluster_inputs).get()
                                   
        # Bring together all the clusters
        star_systems = pd.concat([res_i[0] for res_i in res])
        star_systems.sort_index(inplace=True)
        companions = pd.concat([res_i[1] for res_i in res]) if (self.imf_module.spisea_imf.make_multiples) else None
        # Convert magnitude system if needed
        if (self.evolution_module.photsys_convert is not None) \
                        and (self.evolution_module.bands[0] in props):
            for i,band in enumerate(self.evolution_module.bands):
                star_systems.loc[:,band] += self.evolution_module.photsys_convert[band]      
        if companions is not None:
            companions.sort_index(inplace=True)
            if (self.evolution_module.photsys_convert is not None) \
                        and (self.evolution_module.bands[0] in props):
                for i,band in enumerate(self.evolution_module.bands):
                    companions.loc[:,band] += self.evolution_module.photsys_convert[band]

        star_systems.loc[:,'system_Mass'] = star_systems['Mass']
        if self.imf_module.spisea_imf.make_multiples:
            if len(companions)>0:
                comp_mass_sums = companions.groupby("system_idx")['Mass'].sum()
                primary_idxs = star_systems.index[star_systems['n_companions']>0]
                star_systems.loc[primary_idxs,'system_Mass'] += comp_mass_sums[
                                        primary_idxs]
            companions.reset_index(inplace=True)
        star_systems.reset_index(inplace=True)

        return star_systems, companions

def spisea_props_to_synthpop(tab):
    renames = {'mass_current': 'Mass', 'mass':'iMass', 'logg':'log_g',
                      'N_companions': 'n_companions', 'e':"eccentricity"}
    for col in list(tab.columns):
        if col in renames:
            tab.rename_column(col, renames[col])
    lums = np.array(tab['L'])
    teffs = np.array(tab['Teff'])
    tab['log_L'] = np.log10(lums/const.Lsun_w)
    tab['log_Teff'] = np.log10(teffs)
    tab['log_R'] = np.log10((np.sqrt(lums/(4*np.pi*const.sigma_sb*teffs**4)))/const.Rsun_m)
    nan_mass = np.isnan(tab['Mass'])
    tab['Mass'][nan_mass] = tab['iMass'][nan_mass]
    tab['star_mass'] = tab['Mass']
    if 'systemMass' in tab.columns:
        tab.remove_column('systemMass')
    tab.remove_column('Teff')
    tab.remove_column('L')

    return tab

def generate_spisea_cluster_stars(system_idxs, log_age, mh, feh,
                                  avg_mass_per_star,
                                  block_spisea_prints,
                                  evo_model, atm_func, wd_atm_func,
                                  min_mass, max_mass, bands_obs_str,
                                  imf, ifmr, iso_dir, props,
                                  bands, bbh_frac):
    """
    Separated function to run SPISEA clusters for parallelization
    """
    max_system_idx = -1
    star_systems_list_bin = []
    companions_list_bin = []
    n_bin = len(system_idxs)
    print(f"Starting SPISEA cluster generation for bin log_age={log_age:.2f}"
                        f" [M/H]={mh:.2f} for {n_bin} stars")
    cluster_stars_needed = n_bin
    # Use a minimum mass per cluster of 100.0 so we don't get an error
    generate_mass = np.maximum(cluster_stars_needed*avg_mass_per_star*1.1, 100.0)
    # Loop until we have enough stars
    while cluster_stars_needed > 0:
        with BlockSpiseaPrints(block_prints=block_spisea_prints):
            if evo_model.model_version_name=='COSMIC':
                isochrone = spisea_synthetic.IsochronePhotExternalEvolution(logAge=log_age, AKs=0,
                                    distance=10, metallicity=mh,
                                    evo_model=evo_model, atm_func=atm_func,
                                    wd_atm_func=wd_atm_func, atm_grid_dir=iso_dir,
                                    min_mass=min_mass, max_mass=max_mass,
                                    filters=bands_obs_str)
            else:
                isochrone = spisea_synthetic.IsochronePhot(logAge=log_age, AKs=0,
                                    distance=10, metallicity=mh,
                                    evo_model=evo_model, atm_func=atm_func,
                                    wd_atm_func=wd_atm_func, iso_dir=iso_dir,
                                    min_mass=min_mass, max_mass=max_mass,
                                    filters=bands_obs_str)
            cluster=spisea_synthetic.ResolvedCluster(isochrone, imf, generate_mass,
                                                ifmr=ifmr, keep_low_mass_stars=True)
        star_systems_i = cluster.star_systems
        star_systems_i['system_idx'] = np.arange(len(star_systems_i)) + max_system_idx + 1
        if "companions" in cluster.__dir__():
            companions_i = cluster.companions
            companions_i['system_idx'] += (max_system_idx + 1)
        if len(star_systems_i)>0:
            max_system_idx = star_systems_i['system_idx'].max()
            keep_idx = ((star_systems_i['mass']>min_mass) & (star_systems_i['mass']<max_mass))
            star_systems_i = star_systems_i[keep_idx]
            star_systems_list_bin.append(star_systems_i)
            cluster_stars_needed -= len(star_systems_i)
            if ("companions" in cluster.__dir__()) and (bbh_frac<1.0) \
                        and np.any(star_systems_i['phase']==103):
                bh_indexes = np.where((star_systems_i['phase']==103) & star_systems_i['isMultiple'])[0]
                n_drop = int(np.round(len(bh_indexes)*(1-bbh_frac)))
                bh_drop_indexes = np.random.choice(bh_indexes, size=n_drop, replace=False)
                sys_drop_idxs = star_systems_i['system_idx'][bh_drop_indexes]
                companions_i.remove_rows(np.where(np.isin(companions_i["system_idx"],sys_drop_idxs))[0])
                star_systems_i['N_companions'][bh_drop_indexes] = 0
                star_systems_i['isMultiple'][bh_drop_indexes] = False
                for filt in bands:
                    star_systems_i[filt][bh_drop_indexes] = np.nan
                    
        if "companions" in cluster.__dir__():
            companions_list_bin.append(companions_i)
        else:
            star_systems_i['N_companions'] = 0
            
    star_systems_bin = vstack(star_systems_list_bin)
    if len(companions_list_bin)>0:
        companions_bin = vstack(companions_list_bin)
    else:
        companions_bin = None
    # Drop any excess stars
    if cluster_stars_needed<0:
        star_systems_bin = star_systems_bin[:cluster_stars_needed]
    # Get the data into the expected form
    star_systems_bin = spisea_props_to_synthpop(star_systems_bin)
    
    # Little column adjustments
    star_systems_bin = star_systems_bin[list(props)+['iMass','Mass','system_idx', 'n_companions']]
    star_systems_bin['age'] = 10**log_age / 1e9
    star_systems_bin['Fe/H_initial'] = feh
    # Get companion stars in expected form
    if (companions_bin is not None) and (len(companions_bin)>0):
        # Drop any companions whose systems got dropped
        companions_bin['Fe/H_initial'] = feh
        companions_bin = companions_bin[np.isin(companions_bin['system_idx'], star_systems_bin['system_idx'])]
        companions_bin = spisea_props_to_synthpop(companions_bin)
        companions_bin = companions_bin[list(props)+['iMass','Mass','system_idx', 'eccentricity', 'log_a']]
        if len(companions_bin)>0:
            companions_bin['Fe/H_initial'] = feh
    elif (companions_bin is not None):
        companions_bin = spisea_props_to_synthpop(companions_bin)
        companions_bin = companions_bin[list(props)+['iMass','Mass','system_idx', 'eccentricity', 'log_a']]
        
    # Deal with indexing
    orig_idxs = np.array(star_systems_bin['system_idx'])
    idxs_map = pd.Series(data=system_idxs, index=orig_idxs, dtype=int)
    star_systems_bin['system_idx'] = system_idxs
    if (companions_bin is not None):
        companions_bin['system_idx'] = idxs_map[companions_bin['system_idx']]
        companions_bin = companions_bin.to_pandas(index='system_idx')

    return star_systems_bin.to_pandas(index='system_idx'), companions_bin
