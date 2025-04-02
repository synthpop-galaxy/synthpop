"""
This file contains the Population Class.
Its task is to handle the generation process for each population
according to the dedicated modules for:
1) population density,
2) initial mass function,
3) age distribution,
4) metallicity distribution,
5+6) isochrone systems and interpolator,
7+8) extinction map and law, and
9) kinematics.
"""

__all__ = ["Population"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-07-06"

# Standard Imports
import time
from typing import Tuple, Dict, Union
import pdb
# Non-Standard Imports
import numpy as np
import pandas
from tqdm.auto import tqdm

# Local Imports
# used to allow running as main and importing to another script
try:
    from . import constants as const
except ImportError:
    import constants as const
    import synthpop_utils as sp_utils
    from position import Position
    from star_generator import StarGenerator
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
    from .star_generator import StarGenerator
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


class Population:
    """
    Population subclass for model class.
    Requires a population parameter file to be initialized.

    Attributes
    ----------
    name : string
    l_deg : float
    b_deg : float
    solid_angle : float
    popid : int
    params : dict

    age : age.Age class
    imf : initial_mass_function.InitialMassFunction subclass
    kinematics : kinematics.Kinematics subclass
    metallicity : metallicity.Metallicity subclass
    population_density : population_density.PopulationDensity subclass
    position : position.Position class

    Methods
    -------
    __init__(self,population_file,population_index,isochrones,bands,**kwargs) : None
        initialize the Population class
    assign_subclasses(self, subclasses) : None
        convenience function to assign subsections to Population upon initialization
    update(self,**kwargs) : None
        update population with new kwargs (e.g., new l_deg and b_deg)
    mc_totmass(self,r_inner,r_outer,N) : float [M_sun]
        monte carlo to calculate the total mass in between radius_1 and radius_2 within
        self.solid_angle_sr using N draws
    central_totmass(self,r_inner,r_outer) : float [M_sun]
        use the central density of a slice to estimate the total mass within the slice
    estimate_field(self,) :
        float [M_sun], float [number of stars]
        estimate the total mass and stars in a field
    generate_field(self) : pandas.DataFrame
        generate stars in slices of self.step_size out to self.max_distance
    generate_stars_in_slice(self,r_inner,r_outer,N) :np.array
        generate the initial parameter for N = stars in a slices slice r_inner and r_outer
    check_field(self, TBD...) :
        estimates the number of missing stars in the generated field
    evolve_field(self,data) : float x 4
        estimates the evolved parameter for each star in the field
    """
    def __init__(
            self, pop_params: str, population_index: int, glbl_params: Parameters,
            **positional_kwargs
            ):
        # machinery to import parameter file for a population only given filename
        self.glbl_params = glbl_params
        self.sun = glbl_params.sun
        self.scale_factor = getattr(glbl_params, "scale_factor", 1)

        self.pop_params = pop_params
        self.name = self.pop_params.name
        self.popid = population_index

        logger.create_info_subsection(f"Population {self.popid};  {self.name}")

        logger.info(f"# Initialize Population {self.popid} ({self.name}) from ")
        logger.info(f"pop_file = {self.pop_params._filename!r}")
        json_dump = self.pop_params.model_dump_json(indent=4).replace('\n', '\n    ')
        logger.debug(f"pop_params = {json_dump}")
        
        if hasattr(pop_params, "warp"):
            self.coord_trans = coord_trans.CoordTrans(
                sun=self.glbl_params.sun, **pop_params.warp)
        else:
            self.coord_trans = coord_trans.CoordTrans(
                sun=self.glbl_params.sun, **self.glbl_params.warp)

        # default pointing is Galactic anti-center
        self.l_deg = 180
        self.b_deg = 0
        self.solid_angle_sr = 1e-12
        self.position = Position(
                self.coord_trans, **{**self.glbl_params.window_type, **positional_kwargs}
            )

        # read min and maximum mass from parameters
        self.min_mass = self.glbl_params.mass_lims['min_mass']
        self.max_mass = self.glbl_params.mass_lims['max_mass']
        self.skip_lowmass_stars = self.glbl_params.skip_lowmass_stars

        # get maximum distance and step_size
        self.max_distance = getattr(
                self.pop_params, 'max_distance', self.glbl_params.max_distance)

        self.step_size = getattr(
                self.pop_params, 'distance_step_size', self.glbl_params.distance_step_size)

        self.N_mc_totmass = self.glbl_params.N_mc_totmass

        # star number correction
        self.lost_mass_option = getattr(self.pop_params, "lost_mass_option",
            self.glbl_params.lost_mass_option)

        self.N_av_mass = getattr(self.pop_params, "N_av_mass", self.glbl_params.N_av_mass)

        # add the coordinate transformation module to the kwargs.
        self.pop_params.population_density_kwargs.coord_trans = self.coord_trans
        self.pop_params.kinematics_func_kwargs.coord_trans = self.coord_trans
        self.pop_params.age_func_kwargs.coord_trans = self.coord_trans
        self.pop_params.metallicity_func_kwargs.coord_trans = self.coord_trans
        # add location of sun to kwargs:
        self.pop_params.population_density_kwargs.sun = self.sun
        self.pop_params.kinematics_func_kwargs.sun = self.sun
        self.pop_params.age_func_kwargs.sun = self.sun
        self.pop_params.metallicity_func_kwargs.sun = self.sun

        # placeholder for profiles etc. initialized via the assign_subclasses
        (
            self.population_density, self.imf, self.age, self.metallicity,
            self.kinematics, self.evolution, self.extinction
            ) = self.assign_subclasses()

        # get magnitudes from evolution class
        if isinstance(self.evolution, list):
            self.bands = list(self.evolution[0].bands)
        else:
            self.bands = list(self.evolution.bands)

        # check if main magnitued is in self.bands
        if self.glbl_params.maglim[0] not in self.bands:
            msg = f'{self.glbl_params.maglim[0]}, used as filter for' \
                  f' the magnitude limit is not in {self.bands}'
            logger.critical(msg)
            raise ValueError(msg)

        # check if all bands are in eff_wavelengths:
        if len(not_found := set(self.bands) - self.glbl_params.eff_wavelengths.keys()) != 0:
            raise KeyError(f"Effect Wavelengths for {not_found} are not specified")

        # set wavelength bands and effective wavelength
        self.extinction.set_bands(self.bands, self.glbl_params.eff_wavelengths)

        self.av_mass_corr = None
        self.generator = StarGenerator(
            self.imf, self.age, self.metallicity, self.evolution,
            self.glbl_params, self.position, self.max_mass, logger
            )

    def assign_subclasses(self):
        """initialization of all the subclass
        Returns
        -------
        population_density,
        imf,
        age,
        metallicity,
        kinematics,
        evolution,
        extinction
        """

        evolution = self.get_evolution_class()
        extinction = self.get_extinction_class()
        population_density = self.get_population_density_class()
        imf = self.get_imf_class()
        age = self.get_age_class()
        metallicity = self.get_metallicity_class()
        kinematics = self.get_kinematics_class(population_density)

        return population_density, imf, age, metallicity, kinematics, evolution, extinction

    def get_evolution_class(self):
        evolution_class_config = self.get_evolution_class_config()
        evolutions = []
        # go through  all mass ranges
        # this is to combine different isochrone systems or interpolator
        ev_init = []
        for i, ev_kwargs in enumerate(evolution_class_config):
            # get the Isochrone system

            Isochrone_System = sp_utils.get_subclass(
                EvolutionIsochrones, ev_kwargs,
                initialize=False, population_file=self.pop_params._filename)

            # check if an interpolator is specified
            if hasattr(ev_kwargs, 'interpolator'):
                Interpolator = sp_utils.get_subclass(
                    EvolutionInterpolator, sp_utils.ModuleKwargs(name=ev_kwargs.interpolator),
                    initialize=False, population_file=self.pop_params._filename)
            else:
                Interpolator = None

            # combine Interpolator and Isochrone_System
            Evo = CombineEvolution(Isochrone_System, Interpolator)

            # collect all needed columns and add them to the evolution keywords
            columns = MUST_HAVE_COLUMNS.copy()
            columns.extend(const.REQ_ISO_PROPS)
            columns.extend(self.glbl_params.opt_iso_props)
            # add magnitudes to columns
            if isinstance(self.glbl_params.chosen_bands, dict):
                # 1convert dict into list of tuples
                columns.extend(self.glbl_params.chosen_bands.items())
            else:
                columns.extend(self.glbl_params.chosen_bands)
            ev_kwargs.columns = columns

            ev_init.append(ev_kwargs)
            evolution_i = Evo(
                iso_kwargs=ev_kwargs.init_kwargs,
                int_kwargs=ev_kwargs.init_kwargs,
                logger=logger)

            val = getattr(ev_kwargs, "min_mass", evolution_i.min_mass)
            if val < evolution_i.min_mass:
                msg = f'min mass ({val:.3f}) in the evolution keywords ' \
                      f'is lower than the masses in the Isochrone grid ({evolution_i.min_mass:.3f})'
                logger.critical(msg)
                raise ValueError(msg)

            evolution_i.min_mass = val

            evolution_i.max_mass = getattr(ev_kwargs, "max_mass", evolution_i.max_mass)
            logger.debug(
                "%s : using EvolutionIsochrones subclass '%s' (%.2f-%.2f Msun)", self.name,
                evolution_i.isochrones_name, evolution_i.min_mass, evolution_i.max_mass)
            logger.debug(
                "%s : using EvolutionInterpolator subclass '%s' (%.2f-%.2f Msun)", self.name,
                evolution_i.interpolator_name, evolution_i.min_mass, evolution_i.max_mass)
            evolutions.append(evolution_i)
        logger.log(15, '"evolution_class" : [')
        for ev_kwargs in ev_init:
            msg = f'    {ev_kwargs},'.replace("'", '"')
            msg = msg.replace('False', 'false').replace('True', 'true')
            msg = msg.replace('None', 'null')
            logger.log(15, msg)
        logger.log(15, '   ]')
        if len(evolutions) == 1:
            evolution = evolutions[0]
        else:
            evolution = evolutions
        return evolution

    def get_evolution_class_config(self):
        evolution_class = getattr(self.glbl_params, "evolution_class")
        if evolution_class is None:
            if hasattr(self.pop_params, "evolution_kwargs"):
                logger.debug("read evolution class from Population file")
                evolution_class = self.pop_params.evolution_kwargs
            else:
                msg = "evolution_kwargs is neither defined " \
                      "via the config file nor via the population"
                logger.critical(msg)
                raise ModuleNotFoundError(msg)
        else:
            logger.debug("read evolution class from config file ")
        if not isinstance(evolution_class, list):
            evolution_class = [evolution_class]
        return evolution_class

    def get_extinction_class(self):
        ExtMap = sp_utils.get_subclass(ExtinctionMap, self.glbl_params.extinction_map_kwargs,
                                       initialize=False, population_file=self.pop_params._filename)
        logger.debug(
            [f'initialize Extinction Map with keywords {self.glbl_params.extinction_map_kwargs}',
             ''])
        logger.log(15, f'"extinction_map" :  {self.glbl_params.extinction_map_kwargs}')
        # check if multiple extinction laws are provided
        if not isinstance(self.glbl_params.extinction_law_kwargs, list):
            ExtLaw = sp_utils.get_subclass(ExtinctionLaw, self.glbl_params.extinction_law_kwargs,
                                           initialize=False, population_file=self.pop_params._filename)
        else:
            # collect all extinction laws in a common list
            ExtLaw = [
                sp_utils.get_subclass(ExtinctionLaw, ext_law, initialize=False,
                                      population_file=self.pop_params._filename)
                for ext_law in self.glbl_params.extinction_law_kwargs]
        # combine Extinction and Extinction laws in a combined Class
        Extinction = CombineExtinction(ext_map=ExtMap, ext_law=ExtLaw)
        # initialize extinction
        if not isinstance(self.glbl_params.extinction_law_kwargs, list):
            ext_law_kwargs = [self.glbl_params.extinction_law_kwargs]
        else:
            ext_law_kwargs = self.glbl_params.extinction_law_kwargs
        logger.log(15, '"extinction_law_kwargs" : [')
        for ext_law_kwarg in ext_law_kwargs:
            msg = f'    {ext_law_kwarg},'.replace("'", '"')
            msg = msg.replace('False', 'false').replace('True', 'true')
            msg = msg.replace('None', 'null')
            logger.debug(f'initialize Extinction law with keywords {ext_law_kwarg}')
            logger.log(15, msg)
        logger.log(15, '   ]')
        extinction = Extinction(
            ext_map_kwargs=self.glbl_params.extinction_map_kwargs,
            ext_law_kwargs=self.glbl_params.extinction_law_kwargs,
            logger=logger
        )
        return extinction

    def get_population_density_class(self):
        logger.debug("%s : using PopulationDensity subclass '%s'", self.name,
                     self.pop_params.population_density_kwargs.name)
        population_density = sp_utils.get_subclass(
            PopulationDensity,
            self.pop_params.population_density_kwargs,
            keyword_name='population_density_kwargs',
            population_file=self.pop_params._filename,
        )
        return population_density

    def get_imf_class(self):
        logger.debug("%s : using InitialMassFunction subclass '%s'", self.name,
                     self.pop_params.imf_func_kwargs.name)
        imf = sp_utils.get_subclass(
            InitialMassFunction, self.pop_params.imf_func_kwargs,
            keyword_name='imf_func_kwargs',
            population_file=self.pop_params._filename)
        return imf

    def get_age_class(self):
        logger.debug("%s : using Age subclass '%s'", self.name,
                     self.pop_params.age_func_kwargs.name)
        age = sp_utils.get_subclass(
            Age, self.pop_params.age_func_kwargs,
            keyword_name='age_func_kwargs',
            population_file=self.pop_params._filename)
        return age

    def get_metallicity_class(self):
        logger.debug("%s : using Metallicity subclass '%s'", self.name,
                     self.pop_params.metallicity_func_kwargs.name)
        metallicity = sp_utils.get_subclass(
            Metallicity, self.pop_params.metallicity_func_kwargs,
            keyword_name='metallicity_func_kwargs',
            population_file=self.pop_params._filename)
        return metallicity

    def get_kinematics_class(self, population_density):
        self.pop_params.kinematics_func_kwargs.density_class = population_density
        logger.debug("%s : using Kinematics subclass '%s'", self.name,
                     self.pop_params.kinematics_func_kwargs.name)
        kinematics = sp_utils.get_subclass(
            Kinematics,
            self.pop_params.kinematics_func_kwargs,
            keyword_name='metallicity_func_kwargs',
            population_file=self.pop_params._filename)
        return kinematics

    def set_position(
            self, l_deg: float, b_deg: float, solid_angle: float, solid_angle_unit: str,
            **position_kwargs
            ):
        """
        set a new location and broadcast the new location position and extinction
        Parameters
        ----------
        l_deg : float [degree]
            galactic longitude
        b_deg : float [degree]
            galactic latitude
        solid_angle : float [solid_angle_unit]
            setting size of the cone
        solid_angle_unit : str
            unit for solid_ange
        position_kwargs: dict
            keywords forwarded to the position class
        Returns
        -------

        """

        # set up essential class properties
        self.l_deg = l_deg
        self.b_deg = b_deg
        if solid_angle_unit.lower() in ["sr", "steradian"]:
            self.solid_angle_sr = solid_angle
        elif solid_angle_unit.lower() in ["square degree", "sdeg", "deg^2", "degree^2"]:
            self.solid_angle_sr = solid_angle * (np.pi / 180) ** 2
        else:
            raise ValueError(f"solid angle unit '{solid_angle_unit}' is not supported")

        # flag to mark that coordinates have been set
        self.position.update_location(l_deg, b_deg, self.solid_angle_sr)

        logger.debug(
            f"{self.name} : position is set to "
            f"{self.l_deg: .3f}, {self.b_deg: .3f} with solid angle"
            f"{self.solid_angle_sr:.3e} sr ({solid_angle:.3e} {solid_angle_unit})")

    def mc_totmass(self, r_inner: float, r_outer: float, n_picks: int = 1000) -> float:
        """
        Monte Carlo integration of the total mass in a slice,
        given near and far bounding radii r_inner and r_outer.
        Want to add some catch for large error in the mass,
        then increase value of N.
        Recall that for a Monte Carlo integration:
        Q_N = Volume * (1/N) * sum_N f(x) = V*<f>,
        the expected error in Q_N decreases as 1/sqrt(N).

        Parameters
        ----------
        r_inner : float [kpc]
            inner radius of the slice 
        r_outer : float [kpc]
            outer radius of the slice
        n_picks : int
            number of picks for the  integration

        Returns
        -------
        total_mass : float
            total mass in the slice
        """
        # calculate volume of slice, needs to be kpc^3, r_inner and r_outer in kpc
        volume = (1 / 3) * self.solid_angle_sr * (r_outer ** 3 - r_inner ** 3)

        # MC draws of density
        d, lstar_deg, bstar_deg = self.position.draw_random_point_in_slice(r_inner,
            r_outer, n_picks)[3:]
        r, phi_rad, z = self.coord_trans.dlb_to_rphiz(d, lstar_deg, bstar_deg)
        mean_density = np.mean(self.population_density.density(r, phi_rad, z))
        # then find mass from volume*<density>
        total_mass = volume * mean_density

        return total_mass

    def central_totmass(self, r_inner: float, r_outer: float) -> float:
        """
        Returns the mass using the denisty in the center of a slice

        Parameters
        ----------
        r_inner : float [kpc]
            inner radius of the slice
        r_outer : float [kpc]
            outer radius of the slice

        Returns
        -------
        totmass: float [Msun]
            mass of a slice using the central density
        """

        # calculate volume of slice, needs to be kpc-3, r_inner and r_outer in kpc
        volume = (1 / 3) * self.solid_angle_sr * (r_outer ** 3 - r_inner ** 3)
        r, phi_rad, z = self.coord_trans.dlb_to_rphiz((r_outer + r_inner) / 2, self.l_deg, self.b_deg)
        density = self.population_density.density(r, phi_rad, z)

        totmass = volume * density

        return totmass

    def estimate_field(self, **kwargs) -> Tuple[float, float]:
        """
        estimate the field in the current pointing
        """

        if self.position.l_deg is None:
            msg = ("coordinates where not set."
                   "You must run set_position(l_deg, b_deg, solid_angle_sr)"
                   "before running 'estimate_field' or 'generate_field'")
            logger.critical(msg)
            raise AttributeError(msg)

        average_star_mass = self.imf.average_mass(min_mass=self.min_mass, max_mass=self.max_mass)

        # radii we will step over
        radii = np.arange(0, self.max_distance + self.step_size, self.step_size)
        # find total mass in cone
        # sum over all slices
        total_stellar_mass = sum(self.mc_totmass(inner_radii, outer_radii, n_picks=1000)
            for inner_radii, outer_radii in zip(radii, radii[1:]))
        # find total number of stars
        if self.population_density.density_unit == 'init_mass':
            av_mass_corr = 1
        elif self.lost_mass_option in [1, 2]:
            use_save_av_mass = self.lost_mass_option == 1
            av_mass_corr = self.estimate_average_mass_correction(n_stars=self.N_av_mass,
                use_save_value=use_save_av_mass)
            logger.debug('%s : average mass loss: %.1f%s', self.name, av_mass_corr * 100, '%')
        else:
            if self.pop_params.av_mass_corr is not None:
                av_mass_corr = self.pop_params.av_mass_corr
            elif self.pop_params.n_star_corr is not None:
                av_mass_corr = 1 / self.pop_params.n_star_corr
            else:
                av_mass_corr = 1
        if self.population_density.density_unit == 'number':
            n_star_expected = total_stellar_mass
            total_stellar_mass = n_star_expected * (average_star_mass * av_mass_corr)
        else:
            n_star_expected = total_stellar_mass / (average_star_mass * av_mass_corr)

        if self.scale_factor == 1:
            logger.debug("%s : expected total stellar mass %.3eMsun, expected number of stars %i",
                self.name, total_stellar_mass, n_star_expected)
        else:
            logger.debug("%s : expected total stellar mass %.3eMsun, expected number of stars %i, "
                         "expected scaled number of stars %.f (scale factor %.2f)",
                self.name, total_stellar_mass, n_star_expected,
                n_star_expected / self.scale_factor, self.scale_factor)
        return total_stellar_mass, n_star_expected

    def estimate_average_mass_correction(
            self,
            n_stars: int = 10000,
            use_save_value: bool = False,
            **kwargs
            ) -> float:
        """
        Estimates the ratio between the average evolved mass and initialmass
        Generates and evolve N stars in the cone,
        Evolve them and estimates the average mass

        Parameters
        ----------
        n_stars : int
            numbers of test stars
        use_save_value: bool
            used the stored value if it is already estimated for a previous location

        Returns
        -------
        av_mass_corr : float
            sum(mass_evolved)/sum(mass_initial)
        """

        if use_save_value & (self.av_mass_corr is not None):
            # use previously estimate values
            return self.av_mass_corr

        # generate positions for a test sample
        positions = np.array(
            self.position.draw_random_point_in_slice(0, self.max_distance, n_stars))

        (m_initial, age, met, _, s_props, _, _, not_evolved
        ) = self.generator.generate_star_at_location(positions[0:3].T,
            {const.REQ_ISO_PROPS[0], self.glbl_params.maglim[0]},
            min_mass=self.min_mass, max_mass=self.max_mass)
        # get evolved mass
        mass_evolved = s_props[const.REQ_ISO_PROPS[0]]
        # assume no mass loss for not evolved stars
        mass_evolved[not_evolved] = m_initial[not_evolved]
        # get average mass from imf
        average_m_initial = self.imf.average_mass(min_mass=self.min_mass, max_mass=self.max_mass)
        # estimate mass loss correction
        av_mass_corr = np.mean(mass_evolved) / average_m_initial

        self.av_mass_corr = av_mass_corr

        return av_mass_corr

    def get_mass_loss_for_option(self, lost_mass_option):

        if lost_mass_option == 1:
            # estimate the av_mass_corr in the estimate field
            # and use it for different positions
            av_mass_corr = self.estimate_average_mass_correction(
                n_stars=self.N_av_mass, use_save_value=True
                )

        elif lost_mass_option == 2:
            # estimate the av_mass_corr for each field
            av_mass_corr = self.estimate_average_mass_correction(
                n_stars=self.N_av_mass, use_save_value=False
                )

        elif lost_mass_option == 4:
            # use the av_mass_corr or n_star_corr keyword  specified in the population_config_file
            if self.pop_params.av_mass_corr is not None:
                av_mass_corr = self.pop_params.av_mass_corr
            elif self.pop_params.n_star_corr is not None:
                av_mass_corr = 1 / self.pop_params.n_star_corr
            else:
                av_mass_corr = 0.71

        else:
            av_mass_corr = 1

        return av_mass_corr

    def get_n_star_expected(self, radii, average_imass_per_star, av_mass_corr):
        """ estimates the number of stars in each slice """

        # estimate the mass/numbers in each slice
        mass_per_slice = np.array([self.mc_totmass(radii_inner, radii_outer, self.N_mc_totmass)
            for radii_inner, radii_outer in zip(radii, radii[1:])])
        ################################################################
        #   Translate density into number of generated stars           #
        ################################################################

        if self.population_density.density_unit == "number":
            n_star_expected = mass_per_slice  # density returns n_star_expected
            mass_per_slice = n_star_expected * average_imass_per_star * av_mass_corr

        if self.population_density.density_unit == 'init_mass':
            n_star_expected = mass_per_slice / average_imass_per_star

        else:
            n_star_expected = mass_per_slice / (average_imass_per_star * av_mass_corr)

        return n_star_expected, mass_per_slice

    @staticmethod
    def convert_to_dataframe(
            popid, initial_parameters, m_evolved, final_phase_flag,
            galactic_coordinates, proper_motions, cartesian_coordinates, velocities,
            vr_lsr, extinction_in_map, props, user_props, mags, headers
            ):
        """ convert data to a pandas data_frame"""
        df = pandas.DataFrame(np.column_stack(
            [np.repeat(popid, len(final_phase_flag)),  # pop,
                initial_parameters,  # iMass, age, Fe/H,
                m_evolved, final_phase_flag,  # Mass, In_Final_Phase
                galactic_coordinates, proper_motions,  # Dist, l, b, Vr, mu_l, mu_b,
                cartesian_coordinates, velocities,  # x, y, z,  U, V, W,
                vr_lsr, extinction_in_map,  # VR_LSR, extinction_in_map
                props, user_props, mags  # all specified props and magnitudes
                ]
            ), columns=headers)
        return df

    def generate_field(self) -> Tuple[pandas.DataFrame, Dict]:
        """
        Generate the stars in the field
        estimates the number of stars from a density distribution
        and generates stars based on the individual distributions
        estimates evolved properties by interpolating the isochrones

        Returns
        -------
        population_df : DataFrame
            collected stars for this population
        distribution : dict
            collected distributions, by now only the distance distribution

        """
        if self.position.l_deg is None:
            msg = ("coordinates where not set."
                   "You must run set_position(l_deg, b_deg, solid_angle_sr)"
                   "before running 'estimate_field' or 'generate_field'")
            logger.critical(msg)
            raise AttributeError(msg)

        logger.create_info_subsection(f"Population {self.popid};  {self.name}")

        # array of radii to use
        radii = np.arange(0, self.max_distance + self.step_size, self.step_size)

        # placeholder to collect distributions
        distribution = {}
        # placeholder to collect data frames
        df_list = []

        # collect all the column names
        # required_properties + optional_properties + magnitudes
        headers = const.COL_NAMES + self.glbl_params.col_names + self.bands
        # replace "ExtinctionInMap" with the output of the extinction map
        if 'ExtinctionInMap' in headers:
            extinction_index = headers.index("ExtinctionInMap")
            headers[extinction_index] = self.extinction.A_or_E_type
        # requested properties
        props_list = set(const.REQ_ISO_PROPS + self.glbl_params.opt_iso_props + self.bands)

        # get average_mass from imf
        average_imass_per_star = self.imf.average_mass(min_mass=self.min_mass,
            max_mass=self.max_mass)

        # get initial av_mass_corr # might be estimated on the fly
        av_mass_corr = self.get_mass_loss_for_option(self.lost_mass_option)

        n_star_expected, mass_per_slice = self.get_n_star_expected(
            radii, average_imass_per_star, av_mass_corr)

        if self.lost_mass_option == 3:
            if np.sum(n_star_expected) < self.N_av_mass:
                n_star_expected *= self.N_av_mass / np.sum(n_star_expected)

        mass_limit = np.ones(radii[:-1].shape) * self.min_mass  # new min mass for each
        frac_lowmass = (0., 0.)  # average mass/ fraction of stars

        if (self.glbl_params.maglim[1] > 50) \
                or (not self.skip_lowmass_stars) \
                or (not self.glbl_params.obsmag):
            pass

        elif self.skip_lowmass_stars and isinstance(self.evolution, list):
            logger.warning("skip_lowmass_stars is not implemented for multiple isochrones")

        elif self.skip_lowmass_stars:
            logger.info(f"{self.name} : estimate minimum mass for magnitude limit")
            max_age = self.age.get_maximum_age()
            '''if self.glbl_params.obsmag:
                self.extinction.update_extinction_in_map(radius=radii[:-1])
                _,ext_dict_temp = self.extinction.get_extinctions(self.l_deg*np.ones(len(radii)-1),self.b_deg*np.ones(len(radii)-1),radii[:-1])
                extinction_at_slice_fronts = ext_dict_temp[self.glbl_params.maglim[0]]
            else:
                extinction_at_slice_fronts = None'''
            #print('ext slices',extinction_at_slice_fronts)
            mass_limit = self.evolution.get_mass_min(
                self.glbl_params.maglim[0],
                self.glbl_params.maglim[1], radii[:-1], 
                max_age #, extinction_at_slice_fronts
                )

            mass_limit = np.maximum(mass_limit, self.min_mass)
            mass_limit = np.minimum(mass_limit, self.max_mass)

            frac_lowmass = (
                np.array([self.imf.average_mass(self.min_mass, mm) for mm in mass_limit]),
                ((self.imf.F_imf(mass_limit) - self.imf.F_imf(self.min_mass))
                / (self.imf.F_imf(self.max_mass) - self.imf.F_imf(self.min_mass)))
                )  # average_mass, fraction of stars

        # reduce number of stars by the scale factor and fract_above_min_mass
        total_stars = np.random.poisson(n_star_expected * (1 - frac_lowmass[1]) / self.scale_factor)

        expected_total_imass = n_star_expected * average_imass_per_star

        distribution["distance_distribution"] = np.vstack(
            [(radii[1:] + radii[:-1]) / 2, n_star_expected]).T
        distribution["distance_distribution_comment"] = \
            "pairs of distances in [kpc] and number of stars in the slice"

        logger.info("# From density profile (number density)")

        logger.info(f"expected_total_iMass = {np.sum(expected_total_imass):.4f}")
        logger.info(f"expected_total_eMass = {np.sum(mass_per_slice):.4f}")
        logger.info(f"average_iMass_per_star = {average_imass_per_star:.4f}")
        logger.info(f"mass_loss_correction = {np.sum(av_mass_corr):.4f}")
        logger.info(f"n_expected_stars = {np.sum(n_star_expected):.4f}")
        if self.skip_lowmass_stars:
            logger.info(f"without_lm_stars = {np.sum(n_star_expected * (1 - frac_lowmass[1])):.4f}")
        logger.debug(f"{self.name} : Lost mass option: %s", self.lost_mass_option)
        logger.debug(f"{self.name} : Generate ~{sum(total_stars)} stars")
        self.population_density.average_mass = average_imass_per_star * av_mass_corr
        self.population_density.av_mass_corr = av_mass_corr

        if not self.glbl_params.kinematics_at_the_end:
            logger.info("# Determine velocities when position are generated ")

        ################################################################
        #                  Generate Stars                              #
        ################################################################
        ti = time.time()  # start timer
        missing_stars = total_stars
        loop_counts = 0
        logger.debug("generate stellar properties")
        if self.lost_mass_option==3:
            all_m_initial = []
            all_m_evolved = []
            all_r_inner = []
        opt3_mass_loss_done=False
        use_pbar = np.sum(total_stars)>self.glbl_params.chunk_size
        if use_pbar:
            pbar = tqdm(total=sum(missing_stars))
        while any(missing_stars > 0):
            neg_missing_stars = np.minimum(missing_stars,0)
            missing_stars = np.maximum(missing_stars,0)
            if sum(missing_stars)>self.glbl_params.chunk_size:
                final_expected_loop=False
                idx_cs = np.searchsorted(np.cumsum(missing_stars), self.glbl_params.chunk_size)
                rem_chunk = self.glbl_params.chunk_size - (np.cumsum(missing_stars)[idx_cs-1])*(idx_cs>0)
                missing_stars_chunk = missing_stars * (np.cumsum(missing_stars)<self.glbl_params.chunk_size)
                missing_stars_chunk[idx_cs] = rem_chunk
            else:
                final_expected_loop=True
                missing_stars_chunk = np.copy(missing_stars)

            position, r_inner, proper_motions, velocities, vr_lsr, \
            (m_initial, age, met, ref_mag, s_props, final_phase_flag, 
                inside_grid, not_evolved) = self.generator.generate_stars(radii, 
                missing_stars_chunk, mass_limit, self.do_kinematics, props_list)

            initial_parameters = np.column_stack([m_initial, age, met])

            # extract magnitudes and convert to observed magnitudes if needed
            mags, extinction_in_map = self.extract_magnitudes(
                r_inner, position[:, 3:6], ref_mag, s_props, inside_grid, not_evolved)

            # extract properties
            m_evolved, props, user_props = self.extract_properties(
                m_initial, s_props, const.REQ_ISO_PROPS,
                self.glbl_params.opt_iso_props, inside_grid, not_evolved)

            # Keep track of all stars generated for option 3, until mass loss estimation is complete
            if self.lost_mass_option==3 and not opt3_mass_loss_done:
                all_m_initial += list(m_initial)
                all_m_evolved += list(m_evolved)
                all_r_inner   += list(r_inner)
                if final_expected_loop:
                    missing_stars_evol = self.check_field(
                                radii, average_imass_per_star, np.array(all_m_initial), np.array(all_m_evolved), np.array(all_r_inner),
                                mass_per_slice, frac_lowmass)
                    missing_stars += missing_stars_evol
                    opt3_mass_loss_done=True
            # Subtract out this chunk from the "missing stars"
            missing_stars -= missing_stars_chunk

            # Convert Table to pd.DataFrame
            df = self.convert_to_dataframe(
                self.popid, initial_parameters, m_evolved, final_phase_flag,
                position[:, 3:6], proper_motions, position[:, 0:3], velocities,
                vr_lsr, extinction_in_map, props, user_props, mags, headers
                )

            # add to previous drawn data
            if (self.glbl_params.maglim[-1] != "keep") and (not self.glbl_params.kinematics_at_the_end) and (not self.glbl_params.lost_mass_option==3):
                df = df[df[self.glbl_params.maglim[0]]<self.glbl_params.maglim[1]]
            df_list.append(df)
            loop_counts += 1
            if use_pbar:
                pbar.update(np.sum(missing_stars_chunk))

        # combine the results from the different loops
        if len(df_list) == 0:
            population_df = pandas.DataFrame(columns=headers, dtype=float)
        else:
            population_df = pandas.concat(df_list, ignore_index=True)
        
        # Remove any excess stars
        if self.lost_mass_option==3:
            r_inner=radii[np.searchsorted(radii, population_df['Dist'])-1]
            population_df = self.remove_stars(population_df, r_inner, neg_missing_stars, radii)
            population_df.reset_index(drop=True,inplace=True)

        to = time.time()  # end timer

        ################################################################
        #                       do some logging                        #
        ################################################################

        logger.debug("%s : Generated %i stars in %i cycles (%.2fs)", self.name, len(population_df),
                     loop_counts, to - ti)
        logger.info('# From Generated Field:')

        logger.info(f'generated_stars = {len(population_df)}')
        if len(population_df) != 0:

            logger.info(f'generated_total_iMass = {population_df["iMass"].sum():.4f}')
            gg = population_df.groupby(pandas.cut(population_df.Dist, radii), observed=False)
            if self.skip_lowmass_stars:
                im_incl = (gg["iMass"].sum()
                           + gg.size() * frac_lowmass[0] * frac_lowmass[1]
                           / (1 - frac_lowmass[1])
                           ).sum()

                logger.info(f'generated_total_iMass_incl_lowmass = {im_incl.sum():.4f}')

            logger.info(f'generated_total_eMass = {population_df["Mass"].sum():.4f}')
            if self.skip_lowmass_stars:
                em_incl = (gg["Mass"].sum()
                           + gg.size() * frac_lowmass[0] * frac_lowmass[1]
                           / (1 - frac_lowmass[1])
                           ).sum()
                logger.info(f'generated_total_eMass_incl_lowmass = {em_incl:.4f}')

            logger.info(f'det_mass_loss_corr = '
                        f'{population_df["Mass"].sum() / population_df["iMass"].sum():.4f}')
            if self.skip_lowmass_stars:
                logger.info(f'det_mass_loss_corr_incl_lowmass = {em_incl / im_incl}')

            logger.debug(f'average_mass_per_star = {population_df["Mass"].mean():.4f}')
            if self.skip_lowmass_stars:
                mean_mass = em_incl * ((1 - frac_lowmass[1]) / gg.size()).sum()
                logger.debug(f'average_mass_per_star_incl_lowmass = {mean_mass:.4f}')

        if self.glbl_params.maglim[-1] != 'keep':
            criteria = population_df[self.glbl_params.maglim[0]] < self.glbl_params.maglim[1]
        else:
            criteria = None

        sp_utils.log_basic_statistics(population_df, f"stats_{self.name}", criteria)
        logger.log(25, '# Done')
        logger.flush()

        return population_df, distribution

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
            mass_per_slice: np.ndarray,
            fract_mass_limit: Union[Tuple[float, float], Tuple[np.ndarray, np.ndarray]],
            ) -> np.ndarray:
        """
        Check if stars are missing in a given slice
        Estimates the number of stars that are missing
        Only used for lost_mass_option 3

        Parameters
        ----------
        radii : ndarray
        m_initial : ndarray
        m_evolved : ndarray
        radii_star : ndarray
        mass_per_slice : ndarray
        fract_mass_limit : ndarray or float
        Returns
        -------
        missing_stars : ndarray
            number of missing stars in each slice
        """
        # estimate current initial mass
        m_in = np.array([np.sum(m_initial[radii_star == r]) for r in radii[:-1]])
        # estimate current evolved mass
        m_evo = np.array([np.sum(m_evolved[radii_star == r]) for r in radii[:-1]])

        if self.population_density.density_unit in ['number', 'init_mass']:
            return np.zeros(len(mass_per_slice))

        # option 3.
        # if the sample is large enough otherwise determine the average mass before.
        # otherwise use the determined stars to estimate the average evolved mass:
        # only pick a second time
        average_emass_per_star_mass = (
            fract_mass_limit[0] * fract_mass_limit[1]  # not generated
            + (1 - fract_mass_limit[1]) * np.mean(m_evolved)  # generated
            )
        n_star_expected = mass_per_slice / average_emass_per_star_mass
        total_stars = np.random.poisson(
            n_star_expected * (1 - fract_mass_limit[1]) / self.scale_factor)

        # reduce number of stars by the scale factor
        exist = np.array([np.sum(radii_star == r) for r in radii[:-1]])
        missing_stars = total_stars - exist
        return missing_stars

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
            props, inside_grid=None, not_evolved=None
            ):

        if inside_grid is None:
            inside_grid = np.ones(len(ref_mag), bool)
        if not_evolved is None:
            not_evolved = np.zeros(len(ref_mag), bool)


        mags = np.full((len(ref_mag), len(self.bands)), 9999.)
        extinction_in_map = np.zeros(len(ref_mag))

        dist_module = 5 * np.log10(galactic_coordinates[:, 0] * 100)

        for i, band in enumerate(self.bands):
            mags[:, i] = props[band]

        if self.glbl_params.obsmag:
            ref_mag[inside_grid] += dist_module[inside_grid]
            mags[inside_grid] += dist_module[inside_grid, np.newaxis]

        extinction_in_map, extinction_dict = self.extinction.get_extinctions(
            galactic_coordinates[:, 1],
            galactic_coordinates[:, 2],
            galactic_coordinates[:, 0])

        if self.glbl_params.obsmag:
            ext_mag = extinction_dict.get(self.glbl_params.maglim[0], 0)
            ref_mag[:] += ext_mag
            for i, band in enumerate(self.bands):
                mags[:, i] += extinction_dict.get(band, 0)

        mag_le_limit = ref_mag < self.glbl_params.maglim[1]

        mags[np.logical_not(mag_le_limit)] = np.nan
        mags[np.logical_not(inside_grid)] = np.nan
        mags[not_evolved] = np.nan

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


