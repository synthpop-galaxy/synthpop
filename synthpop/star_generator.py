"""
This file contains the StarGenerator,
It generates stars based on the provided initial distributions,
evolves them and applies the extinction.
"""

__all__ = ["StarGenerator"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2023-03-31"

from typing import Set, Tuple, Dict, Union
import numpy as np
import time
import pandas
import pdb

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
    from modules.evolution import EvolutionIsochrones, EvolutionInterpolator, \
        CombineEvolution, MUST_HAVE_COLUMNS

else:  # continue import when if synthpop is imported
    from . import synthpop_utils as sp_utils
    from .position import Position
    from .synthpop_utils.synthpop_logging import logger
    from .synthpop_utils import coordinates_transformation as coord_trans
    from .synthpop_utils import Parameters
    from .modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from .modules.evolution import EvolutionIsochrones, EvolutionInterpolator, \
        CombineEvolution, MUST_HAVE_COLUMNS


class StarGenerator:
    """
    Star generator object which is used by a Population to generate its member stars.
    
    Attributes
    ----------
    imf_module
    age_module
    met_module
    evolution_module
    kinematics_at_end : bool
        if true, wait until all stars are generated to calculate kinematics
    chunk_size : int
        number of stars to generate per chunk to limit memory use
    ref_band : str
        primary photometric filter for catalog
    position
    max_mass : float
        maximum allowed stellar mass
    """

    def __init__(self, density_module, imf_module, age_module, met_module, evolution_module,
            glbl_params, position, max_mass, ifmr_module, mult_module, bands, logger):
        self.generator_name = 'StarGenerator'
        self.density_module = density_module
        self.imf_module = imf_module
        self.ifmr_module = ifmr_module
        self.mult_module = mult_module
        self.age_module = age_module
        self.met_module = met_module
        if isinstance(evolution_module, list):
            self.evolution_module = evolution_module
        else:
            self.evolution_module = (evolution_module,)
        self.chunk_size = glbl_params.chunk_size
        self.bands = bands
        self.obsmag = glbl_params.obsmag
        self.position=position
        self.max_mass = max_mass
        self.logger = logger
        self.system_mags = False

    def generate_stars(self, radii, missing_stars, mass_limit,
        do_kinematics, props):
        position = np.vstack([
            np.column_stack(self.position.draw_random_point_in_slice(r_inner, r_outer, n_stars, population_density_func=self.density_module.density))
            for r_inner, r_outer, n_stars in zip(radii, radii[1:], missing_stars)
            ])

        min_mass = np.repeat(mass_limit, missing_stars)
        r_inner = np.repeat(radii[:-1], missing_stars)

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

    def generate_star_at_location(self, position, props, min_mass=None, max_mass=None):
        """
        generates stars at the given positions
        """
        n_stars = len(position)
        # Generate base properties: intial mass, age, metallicity
        m_initial = self.imf_module.draw_random_mass(
            min_mass=min_mass, max_mass=max_mass, N=n_stars)
        age = self.age_module.draw_random_age(n_stars)
        met = self.met_module.draw_random_metallicity(
            N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=age)

        # Generate evolved properties
        s_props, final_phase_flag = self.get_evolved_props(m_initial, met, age, props)

        # If assigned, apply IFMR to handle NS and BH evolution
        if self.ifmr_module is not None:
            s_props = self.apply_ifmr(m_initial, met, s_props, final_phase_flag)
                
        # If assigned, generate companions
        if self.mult_module is not None:
            # Adopt metallicity and age of primary; generate init mass and orbits
            pri_ids, m_initial_companions, periods, eccentricities = \
                        self.mult_module.generate_companions(m_initial)
            # Evolve companion stars
            comp_s_props, comp_final_phase_flag = self.get_evolved_props(
                 m_initial_companions, met[pri_ids], age[pri_ids], props)
            # Apply IFMR if present
            if self.ifmr_module is not None:
                comp_s_props = self.apply_ifmr(m_initial_companions,
                    met[pri_ids], comp_s_props, comp_final_phase_flag)

        # Compile star systems table for output
        m_final = s_props['star_mass']
        star_dict = {"iMass": m_initial, "age": age, "Fe/H_initial":met,
                          "n_companions":np.zeros(len(m_initial)),
                          "system_idx": np.arange(len(m_initial)),
                          "Mass": m_final,
                          "system_Mass": m_final}
        star_dict.update(s_props)
        star_systems = pandas.DataFrame.from_dict(star_dict)

        # If assigned, generate companions table and adjust systems table
        if self.mult_module is not None:
            unique_pris, comp_count = np.unique(pri_ids, return_counts=True)
            m_final_companions = comp_s_props['star_mass']
            comp_dict = {"iMass": m_initial_companions, "Mass": m_final_companions,
                         "system_idx": pri_ids, "period": periods,
                         "eccentricity": eccentricities}
            comp_dict.update(comp_s_props)
            companions = pandas.DataFrame.from_dict(comp_dict)
            # Update systems table
            star_systems.loc[unique_pris,"n_companions"] = comp_count
            comp_mass_sums = companions.groupby("system_idx")['Mass'].sum()
            star_systems.loc[comp_mass_sums.index, "system_Mass"] += comp_mass_sums
        else:
            companions = None

        return star_systems, companions

    def get_evolved_props(
            self,
            m_init: np.ndarray,
            met: np.ndarray,
            age: np.ndarray,
            props: Set,
            **kwargs
            ) -> Tuple[np.ndarray, Dict, np.ndarray]:
        """
        evolve the stars using a list of evolution classes given by self.evolution
        each evolution class have a min and max mass range where it should be used.
        the used class are ranked by the order in the list.

        Parameters
        ----------
        m_init : ndarray [Msun]
            list of initial masses for all stars
        met : ndarray [Fe/H]
            list of metallicities for all stars
        age : ndarray [Gyr]
            list of ages for all stars
        props : ndarray
            set of properties which should be interpolated
        kwargs :
            keywords passed to the evolutionary process

        Returns
        -------
        mag: ndarray
            list of the main magnitude
        s_track: dict
            collection of ndarrays for each of the interpolated properties
        inside_grid: ndarray
            used to check if the the star is inside the isochrone grid
        """
        self.logger.debug("Start evolving field")
        ti = time.time()

        # placeholders
        s_track = {p: np.ones(len(m_init)) * np.nan for p in props}
        #mag = np.nan * np.ones(len(m_init))
        inside_grid = np.ones(len(m_init), bool)
        in_final_phase = np.zeros(len(m_init), bool)
        not_performed = np.ones(len(m_init), bool)

        # check if multiple evolution classes are sepecified

        for i, evolution_i in enumerate(self.evolution_module):
            # check if evolution has an atribute which says if numpy arrays can be used
            if hasattr(evolution_i, 'accept_np_arrays'):
                accept_np_arrays = evolution_i.accept_np_arrays
            else:
                accept_np_arrays = True

            # find the stars which fall into the mass range of the current evolution class
            if i != len(self.evolution_module) - 1:
                which = np.where(not_performed & (m_init > evolution_i.min_mass) & (
                        m_init < evolution_i.max_mass))[0]
            else:
                which = np.where(not_performed & (m_init > evolution_i.min_mass))[0]

            # check if there are any stars for this step
            if len(which) == 0:
                continue

            if accept_np_arrays:
                s_props_i, inside_grid_i, in_final_phase_i = evolution_i.get_evolved_props(
                    m_init[which], met[which], age[which], props, **kwargs)

            else:
                # This can be used if self.evolution.get_evolved_props
                # can not handle multiple stars and numpy array:
                m_initial = m_init[which]
                metallicity = met[which]
                age2 = age[which]
                s_props_array, inside_grid_i, in_final_phase_i = np.array([
                    evolution_i.get_evolved_props(
                        m_initial[i], metallicity[i], age2[i] * 1e9, props, **kwargs)
                    for i in range(len(m_initial))
                    ]).T
                # s_props needs to be transformed from an array of dictionaries
                # into a dictionary of numpy arrays
                s_props_i = {
                    key: np.array([i[key] for i in s_props_array])
                    for key in s_props_array[0].keys()}

            # update flags
            inside_grid[which] = inside_grid_i
            in_final_phase[which] = in_final_phase_i

            # update properties
            for key in s_track.keys():
                s_track[key][which] = s_props_i[key]

            # update the list of not performed stars
            not_performed[which] = False
            if not any(not_performed):
                break
        self.logger.debug(f"used time = {time.time() - ti:.2f}s")
        for prop in props:
            if prop == 'star_mass':
                s_track[prop][np.logical_not(inside_grid)] = \
                                m_init[np.logical_not(inside_grid)]
                s_track[prop][not_performed] = m_init[not_performed]
            else:
                s_track[prop][np.logical_not(inside_grid)] = np.nan
                s_track[prop][not_performed] = m_init[not_performed]
        return s_track, in_final_phase

    def apply_ifmr(self, m_init, met, s_props, final_phase_flag):
        """
        Apply the IFMR to catch stars evolved past the grid and make them
        the appropriate dark remnant
        """
        m_compact, m_phase = self.ifmr_module.process_compact_objects(
            m_init[final_phase_flag], met[final_phase_flag])
        #ref_mag[final_phase_flag] = np.nan
        for key in s_props.keys():
            if key=='star_mass':
                s_props[key][final_phase_flag] = m_compact
            elif key=='phase':
                s_props[key][final_phase_flag] = m_phase
            else:
                s_props[key][final_phase_flag] = np.nan
        return s_props




