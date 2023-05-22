
"""
This file contains the StarGenerator,
It generates stars based on the provided initial distributions,
Evolves them and applies the extinction.
"""

__all__ = ["StarGenerator"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2023-03-31"
__license__ = "GPLv3"
__version__ = "0.2.0"

from typing import Set, Tuple, Dict
import numpy as np
import time


class StarGenerator:
    def __init__(
            self, imf_module, age_module, met_module, evolution_module,
            glbl_params, logger):

        self.imf_module = imf_module
        self.age_module = age_module
        self.met_module = met_module
        if isinstance(evolution_module, list):
            self.evolution_module = evolution_module
        else:
            self.evolution_module = (evolution_module,)
        self.chunk_size = glbl_params.chunk_size
        self.ref_band = glbl_params.maglim[0]
        self.logger = logger

    def generate_star_at_location(self, position, props, min_mass=None, max_mass=None):
        """
        generates stars at the given positions
        """
        n_stars = len(position)
        # generate mass
        m_initial = self.imf_module.draw_random_mass(
            min_mass=min_mass, max_mass=max_mass, N=n_stars)

        # generate age
        age = self.age_module.draw_random_age(n_stars)

        # generate  metallicity
        met = self.met_module.draw_random_metallicity(
            N=n_stars, x=position[:, 0], y=position[:, 1], z=position[:, 2], age=age)

        ref_mag, s_props, final_phase_flag, inside_grid, not_evolved = self.get_evolved_props(
            m_initial, met, age, props)

        return m_initial, age, met, ref_mag, s_props, final_phase_flag, inside_grid, not_evolved

    def get_evolved_props(
            self,
            m_init: np.ndarray,
            met: np.ndarray,
            age: np.ndarray,
            props: Set,
            **kwargs
            ) -> Tuple[np.ndarray, Dict, np.ndarray, np.ndarray, np.ndarray]:
        """
        evolve the stars using a list of evolution classes given by 'self.evolution'
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
            collection of numpy ndarrays for each of the interpolated properties
        inside_grid: ndarray
            used to check if the star is inside the isochrone grid
        """
        self.logger.debug("Start evolving field")
        ti = time.time()

        # placeholders
        s_track = {p: np.ones(len(m_init)) * np.nan for p in props}
        mag = np.nan * np.ones(len(m_init))
        inside_grid = np.ones(len(m_init), bool)
        in_final_phase = np.zeros(len(m_init), bool)
        not_performed = np.ones(len(m_init), bool)

        # check if multiple evolution classes are specified

        for i, evolution_i in enumerate(self.evolution_module):
            # check if evolution has an attribute which says if numpy arrays can be used
            if hasattr(evolution_i, 'accept_np_arrays'):
                accept_np_arrays = evolution_i.accept_np_arrays
            else:
                accept_np_arrays = True

            # find the stars which falls into the mass range of the current evolution class
            if i != len(self.evolution_module) - 1:
                which = np.where(not_performed & (m_init > evolution_i.min_mass) & (
                        m_init < evolution_i.max_mass))[0]
            else:
                which = np.where(not_performed & (m_init > evolution_i.min_mass))[0]

            # check if there are any stars for this step
            if len(which) == 0:
                continue
            # loop over bunches of at most chunk_size to reduce memory usage
            if len(which) > self.chunk_size * 6 / 5:
                chunk_size = self.chunk_size
                use_chunks = True
            else:
                chunk_size = len(which) + 1
                use_chunks = False
            count_c = 0

            if use_chunks:
                print(count_c, "/", len(which), end="")

            for which2 in np.array_split(which, len(which) // chunk_size + 1):
                # evolve the stars
                if accept_np_arrays:
                    s_props_i, inside_grid_i, in_final_phase_i = evolution_i.get_evolved_props(
                        m_init[which2], met[which2], age[which2], props, **kwargs)

                else:
                    # This can be used if self.evolution.get_evolved_props
                    # can not handle multiple stars and numpy array:
                    m_initial = m_init[which2]
                    metallicity = met[which2]
                    age2 = age[which2]
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

                # update results to data array
                # update primary magnitude (used for limits etc)
                mag_i = s_props_i.get(self.ref_band)
                mag[which2] = mag_i

                # update flags
                inside_grid[which2] = inside_grid_i
                in_final_phase[which2] = in_final_phase_i

                # update properties
                for key in s_track.keys():
                    s_track[key][which2] = s_props_i[key]

                count_c += len(which2)
                if use_chunks:
                    print("\r", count_c, "/", len(which), end='')

            # update the list of not performed stars
            not_performed[which] = False
            if use_chunks:
                print('')
            # check if anything left to do
            if not any(not_performed):
                break
        self.logger.debug(f"used time = {time.time() - ti:.2f}s")

        return mag, s_track, in_final_phase, inside_grid, not_performed
