
"""
This file contains the in-progress Spisea alternative to StarGenerator,
It generates stars based on the provided initial distributions,
evolves them and applies the extinction.
"""

__all__ = ["SpiseaGenerator"]
__author__ = "M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2023-03-31"
__license__ = "GPLv3"
__version__ = "0.2.0"

from typing import Set, Tuple, Dict
import numpy as np
import time

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
    from modules.age import Age
    from modules.initial_mass_function import InitialMassFunction
    from modules.kinematics import Kinematics
    from modules.metallicity import Metallicity
    from modules.population_density import PopulationDensity

else:  # continue import when if synthpop is imported
    from . import synthpop_utils as sp_utils
    from .position import Position
    from .synthpop_utils.synthpop_logging import logger
    from .synthpop_utils import coordinates_transformation as coord_trans
    from .synthpop_utils import Parameters
    from .modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from .modules.evolution import EvolutionIsochrones, EvolutionInterpolator, \
        CombineEvolution, MUST_HAVE_COLUMNS
    from .modules.age import Age
    from .modules.initial_mass_function import InitialMassFunction
    from .modules.kinematics import Kinematics
    from .modules.metallicity import Metallicity
    from .modules.population_density import PopulationDensity


class SpiseaGenerator:
    def __init__(self, imf_module, age_module, met_module, evolution_module,
            glbl_params, position, max_mass, logger):

        self.imf_module = imf_module
        self.age_module = age_module
        self.met_module = met_module
        if isinstance(evolution_module, list):
            self.evolution_module = evolution_module
        else:
            self.evolution_module = (evolution_module,)
        self.kinematics_at_the_end = glbl_params.kinematics_at_the_end
        self.chunk_size = glbl_params.chunk_size
        self.ref_band = glbl_params.maglim[0]
        self.position=position
        self.max_mass = max_mass
        self.logger = logger

    def generate_stars(self, radii, missing_stars, mass_in_slice, mass_limit, 
        do_kinematics, props):
        position = np.vstack([
            np.column_stack(self.position.draw_random_point_in_slice(r_inner, r_outer, n_stars))
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
        return position, r_inner, proper_motions, velocities, vr_lsr, \
            self.generate_star_at_location(
            position[:, 0:3], props, min_mass, self.max_mass)

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
            N=n_stars, x=position[:,0], y=position[:,1], z=position[:,2], age=age)

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
            collection of ndarratys for each of the interpolated properties
        inside_grid: ndarray
            used to check if the the star is inside the isochrone grid
        """
        self.logger.debug("Start evolving field")
        ti = time.time()

        # placeholders
        s_track = {p: np.ones(len(m_init)) * np.nan for p in props}
        mag = np.nan * np.ones(len(m_init))
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

            # find the stars which falles into the mass range of the current evolution class
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
            if use_chunks: print('')
            # check if anything left to do
            if not any(not_performed):
                break
        self.logger.debug(f"used time = {time.time() - ti:.2f}s")

        return mag, s_track, in_final_phase, inside_grid, not_performed

    @staticmethod
    def remove_stars(
            df: pandas.DataFrame,
            radii_star: np.ndarray,
            missing_stars: np.ndarray,
            radii: np.ndarray
            ) -> pandas.DataFrame:
        """
        Removes stars form data frame in the corresponding slice
        if missing_stars is < 0
        """
        preserve = np.ones(len(radii_star), bool)
        for r, n in zip(radii, missing_stars):
            if n < 0:
                t = preserve[radii_star == r]
                t[n:] = False
                preserve[radii_star == r] = t
        return df[preserve]

    def check_field(
            self,
            radii: np.ndarray,
            average_imass_per_star: float,
            m_initial: np.ndarray,
            m_evolved: np.ndarray,
            radii_star: np.ndarray,
            loop_counts: int,
            mass_per_slice: np.ndarray,
            fract_mass_limit: Union[Tuple[float, float], Tuple[np.ndarray, np.ndarray]],
            ) -> np.ndarray:
        """
        Check if stars are missing in a given slice
        Estimates the number of stars that are missing

        Parameters
        ----------
        radii : ndarray
        m_initial : ndarray
        m_evolved : ndarray
        radii_star : ndarray
        loop_counts : int
        mass_per_slice : ndarray
        fract_mass_limit : ndarray or float
        Returns
        -------
        missing_stars : ndarray
            number of missing stars in each slice
        """

        # estimate current initial mass
        m_in = np.array([np.sum(m_initial[radii_star == r]) for r in radii[:-1]])
        # estimate current evoled mass
        m_evo = np.array([np.sum(m_evolved[radii_star == r]) for r in radii[:-1]])

        if self.population_density.density_unit in ['number', 'init_mass']:
            return np.zeros(1)

        # Option 1, 2 and 6 does not make use of check_field
        elif self.lost_mass_option == 3:
            # option 3.
            # if the sample is large enough otherwise determine the average mass before.
            # otherwise use the determined stars to estimate the average evolved mass:
            # only pick a second time
            if loop_counts > 0:
                return np.zeros(1)

            average_emass_per_star_mass = (
                fract_mass_limit[0] * fract_mass_limit[1]  # not generated
                + (1 - fract_mass_limit) * np.mean(m_evolved)  # generated
                )
            n_star_expected = mass_per_slice / average_emass_per_star_mass
            total_stars = np.random.poisson(
                n_star_expected * (1 - fract_mass_limit[1]) / self.scale_factor)

            # reduce number of stars by the scale factor
            exist = np.array([np.sum(radii_star == r) for r in radii[:-1]])
            missing_stars = total_stars - exist
            return missing_stars

        else:
            return np.zeros(1)

    def do_kinematics(
            self, dist: np.ndarray, star_l_deg: np.ndarray, star_b_deg: np.ndarray,
            x: np.ndarray, y: np.ndarray, z: np.ndarray, **kwargs
            ) -> Tuple[np.ndarray, ...]:
        """
        calls the generation of velocities

        Parameters
        ----------
        dist, star_l_deg, star_b_deg: ndarray
            Spherical Coordinates, used to transform from cartesian
            velocities to galactic velocities
        x, y, z: ndarray
            Cartesian coordinates, used to generate random Veloities
        kwargs : dict
            Keyword arguments passed to draw_random_velocity.

        Returns
        -------
        u, v, w : ndarray
            cartesian velocities
        vr: ndarray
            radial velocity respect to the center of the galaxie
        mu_l, mu_b : ndarray
            galactic proper motion
        rv : ndarray
            radial velocity respect to the earth
        """
        u: ndarray
        v: ndarray
        w: ndarray
        # draw random velocities,
        u, v, w = self.kinematics.draw_random_velocity(x, y, z,
            density_class=self.population_density, **kwargs)

        # translate u, v, w into mu_l and mu_b_vr
        vr, mu_l, mu_b = self.coord_trans.uvw_to_vrmulb(star_l_deg, star_b_deg, dist, u, v, w)

        # correct for motion of the sun
        # following Beaulieu et al. (2000)
        vr_lsr = (
            vr
            + self.sun.v_lsr * np.sin(star_l_deg * np.pi / 180)
            * np.sin(star_b_deg * np.pi / 180)
            + self.sun.v_pec * (
                np.sin(star_b_deg * np.pi / 180) * np.sin(self.sun.b_apex_deg * np.pi / 180)
                + np.cos(star_b_deg * np.pi / 180) * np.cos(self.sun.b_apex_deg * np.pi / 180)
                * np.cos(star_l_deg * np.pi / 180 - self.sun.l_apex_deg * np.pi / 180)
                )
            )

        return u, v, w, vr, mu_l, mu_b, vr_lsr

    def extract_magnitudes(
            self, radii_inner, galactic_coordinates, ref_mag,
            props, inside_grid=None
            ):
        if inside_grid is None:
            inside_grid = np.ones(len(ref_mag), bool)

        mags = np.full((len(ref_mag), len(self.bands)), 9999.)
        extinction_in_map = np.zeros(len(ref_mag))

        dist_module = 5 * np.log10(galactic_coordinates[:, 0] * 100)

        for i, band in enumerate(self.bands):
            mags[:, i] = props[band]

        if self.glbl_params.obsmag:
            ref_mag[inside_grid] += dist_module[inside_grid]
            mags[inside_grid] += dist_module[inside_grid, np.newaxis]

        for ri in np.unique(radii_inner):
            current_slice = radii_inner == ri

            self.extinction.update_extinction_in_map(radius=ri)
            extinction_in_map[current_slice], extinction_dict = self.extinction.get_extinctions(
                galactic_coordinates[current_slice, 1],
                galactic_coordinates[current_slice, 2],
                galactic_coordinates[current_slice, 0])

            if self.glbl_params.obsmag:
                ext_mag = extinction_dict.get(self.glbl_params.maglim[0], 0)
                ref_mag[current_slice] += ext_mag
                for i, band in enumerate(self.bands):
                    mags[current_slice, i] += extinction_dict.get(band, 0)

        mag_le_limit = ref_mag < self.glbl_params.maglim[1]
        mags[np.logical_not(mag_le_limit)] = np.nan

        return mags, extinction_in_map

    @staticmethod
    def extract_properties(
            m_initial, props, req_keys, user_keys,
            inside_grid=None, not_evolved=None
            ):
        # default masks
        if inside_grid is None:
            inside_grid = np.ones(len(ref_mag), bool)
        if not_evolved is None:
            not_evolved = np.zeros(len(ref_mag), bool)

        # setup numpy array to store data
        default_props = np.zeros((len(m_initial), len(req_keys)))
        user_props = np.zeros((len(m_initial), len(user_keys)))

        # store key into default_props and user_props
        for i, key in enumerate(req_keys):
            default_props[:, i] = props[key]
        for i, key in enumerate(user_keys):
            user_props[:, i] = props[key]

        # replace data outside grid with nan
        default_props[np.logical_not(inside_grid), 1:] = np.nan
        default_props[not_evolved, 1:] = np.nan
        default_props[not_evolved, 0] = m_initial[not_evolved]

        user_props[np.logical_not(inside_grid)] = np.nan
        user_props[not_evolved] = np.nan

        return default_props[:, 0], default_props[:, 1:], user_props

