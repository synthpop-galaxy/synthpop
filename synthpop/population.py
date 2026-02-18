"""
This file contains the Population Class.
Its task is to handle the generation process for each population
according to the dedicated modules for:
1) population density,
2) initial mass function,
3) age distribution,
4) metallicity distribution,
5+6) isochrone systems and interpolator,
7) initial-final mass relation,
8) multiplicity,
9+10) extinction map and law, and
11) kinematics.
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
import pandas as pd
from tqdm.auto import tqdm
import warnings

# Local Imports
# used to allow running as main and importing to another script
try:
    from . import constants as const
    from . import synthpop_utils as sp_utils
    from .star_generator import StarGenerator
    from .synthpop_utils.synthpop_logging import logger
    from .synthpop_utils.utils_functions import combine_system_mags, get_primary_mags
    from .synthpop_utils import coordinates_transformation as coord_trans
    from .synthpop_utils import Parameters
    from .modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from .modules.evolution import EvolutionIsochrones, EvolutionInterpolator, \
        CombineEvolution, MUST_HAVE_COLUMNS
    from .modules.age import Age
    from .modules.initial_mass_function import InitialMassFunction
    from .modules.initial_final_mass_relation import InitialFinalMassRelation
    from .modules.kinematics import Kinematics
    from .modules.metallicity import Metallicity
    from .modules.population_density import PopulationDensity
    from .modules.multiplicity import Multiplicity
except ImportError:
    import constants as const
    import synthpop_utils as sp_utils
    from star_generator import StarGenerator
    from synthpop_utils.synthpop_logging import logger
    from synthpop_utils.utils_functions import combine_system_mags, get_primary_mags
    from synthpop_utils import coordinates_transformation as coord_trans
    from synthpop_utils import Parameters
    from modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from modules.evolution import EvolutionIsochrones, EvolutionInterpolator, \
        CombineEvolution, MUST_HAVE_COLUMNS
    from modules.age import Age
    from modules.initial_mass_function import InitialMassFunction
    from modules.initial_final_mass_relation import InitialFinalMassRelation
    from modules.kinematics import Kinematics
    from modules.metallicity import Metallicity
    from modules.population_density import PopulationDensity
    from modules.multiplicity import Multiplicity

class Population:
    """
    Population subclass for model class.
    Requires a population parameter file to be initialized.

    Attributes
    ----------
    name : string
    l_deg : float
    b_deg : float
    popid : int
    params : dict

    age : age.Age class
    imf : initial_mass_function.InitialMassFunction subclass
    kinematics : kinematics.Kinematics subclass
    metallicity : metallicity.Metallicity subclass
    population_density : population_density.PopulationDensity subclass
    position : position.Position class
    ifmr : initial_final_mass_ratio.InitialFinalMassRatio subclass or None
    mult : multiplicity.Multiplicity subclass or None

    Methods
    -------
    __init__(self,population_file,population_index,isochrones,bands,**kwargs) : None
        initialize the Population class
    assign_subclasses(self, subclasses) : None
        convenience function to assign subsections to Population upon initialization
    update(self,**kwargs) : None
        update population with new kwargs (e.g., new l_deg and b_deg)
    estimate_field(self,) :
        float [M_sun], float [number of stars]
        estimate the total mass and stars in a field
    generate_field(self) : pd.DataFrame
        generate stars out to self.max_distance
    check_field(self, TBD...) :
        estimates the number of missing stars in the generated field
    evolve_field(self,data) : float x 4
        estimates the evolved parameter for each star in the field
    """
    def __init__(
            self, pop_params: str, population_index: int, glbl_params: Parameters):
        # machinery to import parameter file for a population only given filename
        self.glbl_params = glbl_params
        self.sun = glbl_params.sun
        self.scale_factor = getattr(glbl_params, "scale_factor", 1)

        self.pop_params = pop_params
        self.name = self.pop_params.name
        self.popid = population_index
        self.obsmag = glbl_params.obsmag

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

        # read min and maximum mass from parameters
        self.min_mass = self.glbl_params.mass_lims['min_mass']
        self.max_mass = self.glbl_params.mass_lims['max_mass']
        self.skip_lowmass_stars = self.glbl_params.skip_lowmass_stars

        # get maximum distance and step_size
        self.max_distance = getattr(
                self.pop_params, 'max_distance', self.glbl_params.max_distance)

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
            self.kinematics, self.evolution, self.extinction, self.ifmr, self.mult
            ) = self.assign_subclasses()

        # get magnitudes from evolution class
        if isinstance(self.evolution, list):
            self.bands = list(self.evolution[0].bands)
            self.eff_wavelengths = dict(self.evolution[0].eff_wavelengths)
        else:
            self.bands = list(self.evolution.bands)
            self.eff_wavelengths = dict(self.evolution.eff_wavelengths)

        # check if main magnitued is in self.bands
        if (self.glbl_params.maglim is not None) and (self.glbl_params.maglim[0] not in self.bands):
            msg = f'{self.glbl_params.maglim[0]}, used as filter for' \
                  f' the magnitude limit is not in {self.bands}'
            logger.critical(msg)
            raise ValueError(msg)

        # check if all bands are in eff_wavelengths:
        if len(not_found := set(self.bands) - self.eff_wavelengths.keys()) != 0:
            raise KeyError(f"Effect Wavelengths for {not_found} are not specified")

        # set wavelength bands and effective wavelength
        if self.extinction is not None:
            self.extinction.set_bands(self.bands, self.eff_wavelengths)
            self.extinction.validate_extinction()

        self.av_mass_corr = None
        if self.glbl_params.star_generator=="SpiseaGenerator":
            if self.skip_lowmass_stars:
                logger.warning("Setting skip_lowmass_stars to False for SpiseaGenerator.")
                self.skip_lowmass_stars = False
            try:
                from spisea_generator import SpiseaGenerator
            except:
                from .spisea_generator import SpiseaGenerator
            self.generator = SpiseaGenerator(self.population_density,
                self.imf, self.age, self.metallicity, self.evolution,
                self.glbl_params, self.max_mass, 
                self.ifmr, self.mult, self.bands, logger
                )
        else:
            self.generator = StarGenerator(self.population_density,
                self.imf, self.age, self.metallicity, self.evolution,
                self.glbl_params, self.max_mass, 
                self.ifmr, self.mult, self.bands, logger
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
        extinction,
        ifmr,
        multiplicity
        """

        evolution = self.get_evolution_class()
        extinction = self.get_extinction_class()
        population_density = self.get_population_density_class()
        imf = self.get_imf_class()
        age = self.get_age_class()
        metallicity = self.get_metallicity_class()
        kinematics = self.get_kinematics_class(population_density)
        ifmr = self.get_ifmr_class()
        mult = self.get_mult_class()

        return population_density, imf, age, metallicity, kinematics, evolution, extinction, ifmr, mult

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
        if (ExtLaw is None) or (ExtMap is None):
            logger.warning('Extinction Map and/or Law is set to None -- no extinction will be applied.')
            return None
        else:
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

    def get_ifmr_class(self):
        ifmr = sp_utils.get_subclass(InitialFinalMassRelation, self.glbl_params.ifmr_kwargs,
                        initialize=True, population_file=self.pop_params._filename)
        return ifmr

    def get_mult_class(self):
        mult = sp_utils.get_subclass(Multiplicity, self.glbl_params.multiplicity_kwargs,
                        initialize=True, population_file=self.pop_params._filename)
        return mult

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
            self, l_deg: float, b_deg: float, field_shape: str, field_scale: float or tuple or np.ndarray, 
            field_scale_unit: str
            ):
        """
        set a new location and broadcast the new location position and extinction
        Parameters
        ----------
        l_deg : float [degree]
            galactic longitude
        b_deg : float [degree]
            galactic latitude
        field_shape : str
            shape of field on-sky: circle or box
        field_scale : float 
            setting size of the cone or box
        field_scale_unit : str
            unit for field_scale
        Returns
        -------

        """

        # set up essential class properties
        self.l_deg = l_deg
        self.b_deg = b_deg
        self.field_shape=field_shape
        if hasattr(field_scale, '__len__'):
            field_scale = np.array(field_scale)
        if field_scale_unit.lower() in ["rad", "radian", 'radians']:
            self.field_scale_deg = field_scale / (np.pi / 180)
        elif field_scale_unit.lower() in ["deg", "degree", 'degrees']:
            self.field_scale_deg = field_scale
        elif field_scale_unit.lower() in ['am', 'arcminute', 'arcmin', 'arcminutes']:
            self.field_scale_deg = field_scale / 60
        elif field_scale_unit.lower() in ['as', 'arcsecond', 'arcsec', 'arcseconds']:
            self.field_scale_deg = field_scale / 60**2
        else:
            raise ValueError(f"field_scale unit '{field_scale_unit}' is not supported")

        # flag to mark that coordinates have been set
        self.population_density.update_location(l_deg, b_deg, field_shape, self.field_scale_deg,
                                                self.max_distance)

        logger.debug(f"{self.name} : position is set to {l_deg: .3f}, {b_deg: .3f}")
        if field_shape=='circle':
            logger.debug(f"{self.field_scale_deg} deg radius circle. ({field_scale} {field_scale_unit})")
        if field_shape=='box':
            logger.debug(f"{self.population_density.l_length_deg}, {self.population_density.b_length_deg} degree l, b length box")

        if self.extinction is not None:
            assert not np.any(np.isnan(self.extinction.get_extinctions(np.array([l_deg]), np.array([b_deg]), np.array([self.max_distance]))[0])), \
                fr"{self.extinction.extinction_map_name} not valid in direction l_deg={l_deg}, b_deg={b_deg} " \
                f"at max distance {self.max_distance}. Check the map's sky coverage."

    def estimate_field(self, **kwargs) -> Tuple[float, float]:
        """
        Estimate the total mass and number of stars in the current field
        """

        if self.population_density.l_deg is None:
            msg = ("coordinates where not set."
                   "You must run set_position(l_deg, b_deg, field_shape, field_scale, field_scale_unit)"
                   "before running 'estimate_field' or 'generate_field'")
            logger.critical(msg)
            raise AttributeError(msg)

        average_star_mass = self.imf.average_mass(min_mass=self.min_mass, max_mass=self.max_mass)

        # total mass is computed in population_density.update_location
        total_stellar_mass = self.population_density.total_mass
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
        Estimates the ratio between the average evolved mass and initial mass by
        generating and evolving n_stars stars in the field andevolving them.

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

        logger.info(f"Evolving test set from population {self.popid} to estimate average initial->final mass ratio")
        # generate positions for a test sample
        positions = np.array(
            self.population_density.draw_random_positions(n_stars))

        star_sample, _ = self.generator.generate_star_at_location(positions[0:4],
            {'star_mass'},
            min_mass=self.min_mass, max_mass=self.max_mass)
        # get evolved mass
        average_m_evolved = star_sample['system_Mass'].mean()
        # get average mass
        average_m_initial = star_sample['iMass'].mean()
        # estimate mass loss correction
        av_mass_corr = average_m_evolved / average_m_initial

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

    def get_n_star_expected(self, average_imass_per_star, av_mass_corr):
        """ estimates the number of stars needed """

        if self.population_density.density_unit == "number":
            n_star_expected = self.population_density.total_mass  # density returns n_star_expected

        if self.population_density.density_unit == 'init_mass':
            n_star_expected = self.population_density.total_mass / average_imass_per_star

        else:
            n_star_expected = self.population_density.total_mass / (average_imass_per_star * av_mass_corr)

        return n_star_expected

    # @staticmethod
    # def convert_to_dataframe(
    #         popid, initial_parameters, m_evolved, final_phase_flag,
    #         galactic_coordinates, proper_motions, cartesian_coordinates, velocities,
    #         vr_lsr, extinction_in_map, props, user_props, mags, headers
    #         ):
    #     """ convert data to a pandas data_frame"""
    #     df = pd.DataFrame(np.column_stack(
    #         [np.repeat(popid, len(final_phase_flag)),  # pop,
    #             initial_parameters,  # iMass, age, Fe/H,
    #             m_evolved, final_phase_flag,  # Mass, In_Final_Phase
    #             galactic_coordinates, proper_motions,  # Dist, l, b, Vr, mu_l, mu_b,
    #             cartesian_coordinates, velocities,  # x, y, z,  U, V, W,
    #             vr_lsr, extinction_in_map,  # VR_LSR, extinction_in_map
    #             props, user_props, mags  # all specified props and magnitudes
    #             ]
    #         ), columns=headers)
    #     return df

    def generate_field(self) -> pd.DataFrame:
        """
        Generate the stars in the field. The number of stars is determined by
        the density distribution, and the individual evolved properties come
        from isochrone interpolation. 

        Returns
        -------
        population_df : DataFrame
            collected stars for this population
        population_companions_df : DataFrame
            companions for the population_df stars (None if no multiplicity)

        """
        if self.population_density.l_deg is None:
            msg = ("Coordinates were not set."
                   "You must run set_position(l_deg, b_deg, field_shape, field_scale, field_scale_unit)"
                   "before running 'estimate_field' or 'generate_field'")
            logger.critical(msg)
            raise AttributeError(msg)

        logger.create_info_subsection(f"Population {self.popid};  {self.name}")

        # array of radii to use
        radii = np.linspace(0, self.max_distance, int(round(self.max_distance/0.1))+1)

        # placeholder to collect data frames
        df_list = []
        comp_df_list = []

        # collect all the column names
        # required_properties + optional_properties + magnitudes
        headers = const.REQ_COL_NAMES + self.glbl_params.opt_iso_props + self.bands
        # replace "ExtinctionInMap" with the output of the extinction map
        if self.extinction is not None:
            headers.append(self.extinction.A_or_E_type)
        # requested properties
        props_list = set(const.REQ_ISO_PROPS + self.glbl_params.opt_iso_props + self.bands)

        # get average_mass from imf
        average_imass_per_star = self.imf.average_mass(min_mass=self.min_mass,
            max_mass=self.max_mass)

        # get initial av_mass_corr # might be estimated on the fly
        av_mass_corr = self.get_mass_loss_for_option(self.lost_mass_option)

        n_star_expected = self.get_n_star_expected(average_imass_per_star, av_mass_corr)

        if (self.lost_mass_option == 3) and (self.population_density.density_unit != 'number') and (n_star_expected>0):
            if n_star_expected < self.N_av_mass:
                n_star_expected = self.N_av_mass

        mass_limit = np.ones(radii[:-1].shape) * self.min_mass  # new min mass for each
        frac_lowmass = (0., 0.)  # average mass/ fraction of stars

        if (self.glbl_params.maglim is None) \
                or (self.glbl_params.maglim[1] > 50) \
                or (not self.skip_lowmass_stars) \
                or (not self.glbl_params.obsmag):
            pass

        elif self.skip_lowmass_stars and isinstance(self.evolution, list):
            logger.warning("skip_lowmass_stars is not implemented for multiple isochrones")
            self.skip_lowmass_stars = False

        elif self.skip_lowmass_stars and self.mult is not None:
            logger.warning("skip_lowmass_stars is not implemented for models with stellar multiplicity")
            self.skip_lowmass_stars = False

        elif self.skip_lowmass_stars:
            raise NotImplementedError("SORRY I BROKE THIS OPTION TEMPORARILY - skip_lowmass_stars not available")
            logger.info(f"{self.name} : estimate minimum mass for magnitude limit")
            max_age = self.age.get_maximum_age()
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
            #pdb.set_trace()
        # reduce number of stars by the scale factor and fract_above_min_mass
        total_stars = np.random.poisson(n_star_expected / self.scale_factor)

        expected_total_imass = n_star_expected * average_imass_per_star

        logger.info("# From density profile (number density)")

        logger.info(f"expected_total_iMass = {expected_total_imass:.4f}")
        logger.info(f"expected_total_eMass = {self.population_density.total_mass:.4f}")
        logger.info(f"average_iMass_per_star = {average_imass_per_star:.4f}")
        logger.info(f"mass_loss_correction = {np.sum(av_mass_corr):.4f}")
        logger.info(f"n_expected_stars = {n_star_expected:.4f}")
        # no longer meaningful when slicing is ONLY used for mass cut, not count
        # if self.skip_lowmass_stars:
        #     logger.info(f"without_lm_stars = {np.sum(n_star_expected * (1 - frac_lowmass[1])):.4f}")
        logger.debug(f"{self.name} : Lost mass option: %s", self.lost_mass_option)
        logger.debug(f"{self.name} : Generate ~{total_stars} stars")
        self.population_density.average_mass = average_imass_per_star * av_mass_corr
        self.population_density.av_mass_corr = av_mass_corr

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
        use_pbar = (np.sum(total_stars)>self.glbl_params.chunk_size) * (self.generator.generator_name!='SpiseaGenerator')
        if use_pbar:
            pbar = tqdm(total=missing_stars)
        neg_missing_stars = np.minimum(missing_stars,0)
        gen_missing_stars = np.maximum(missing_stars,0)
        current_max_id = -1
        while (missing_stars > 0):
            if (gen_missing_stars>self.glbl_params.chunk_size) and not (self.generator.generator_name=='SpiseaGenerator'):
                final_expected_loop=False
                gen_stars_chunk = self.glbl_params.chunk_size
            else:
                final_expected_loop=True
                gen_stars_chunk = gen_missing_stars

            df, comp_df = self.generate_stars(gen_stars_chunk, 
                mass_limit, props_list, radii=radii)

            # Make sure main df has selection option for primary or system mags
            if self.extinction is not None:
                df, comp_df = self.apply_extinction(df, comp_df)
            if self.glbl_params.combine_system_mags and (comp_df is not None) \
                        and (not self.generator.system_mags):
                df = combine_system_mags(df, comp_df, self.bands)
            if (not self.glbl_params.combine_system_mags) and (comp_df is not None) \
                        and (self.generator.system_mags):
                df = get_primary_mags(df, comp_df, self.bands)
                
            # Keep track of all stars generated for option 3, until mass loss estimation is complete
            if self.lost_mass_option==3 and not opt3_mass_loss_done:
                all_m_initial += list(star_set['iMass'])
                all_m_evolved += list(star_set['system_Mass'])
                #all_r_inner   += list(star_set['r_inner'])
                if final_expected_loop:
                    missing_stars = self.check_field(
                                average_imass_per_star, np.array(all_m_initial), np.array(all_m_evolved),
                                self.population_density.total_mass, frac_lowmass)
                    opt3_mass_loss_done=True
                else:
                    missing_stars -= gen_stars_chunk
            # Subtract out this chunk from the "missing stars"
            else:
                missing_stars -= gen_stars_chunk

            # add to previous drawn data
            if (self.glbl_params.maglim is not None) and (not self.glbl_params.lost_mass_option==3):
                df = df[df[self.glbl_params.maglim[0]]<self.glbl_params.maglim[1]]
            if (len(df)>0) and (self.mult is not None):
                comp_df = comp_df[np.isin(comp_df['system_idx'].to_numpy(), df['system_idx'].to_numpy())]
            elif (len(df)==0) and (self.mult is not None):
                comp_df.drop(comp_df.index, inplace=True)
                
            df.loc[:,'system_idx'] += (current_max_id + 1)
            if self.mult is not None:
                comp_df.loc[:,'system_idx'] += (current_max_id + 1)
            if len(df)>0:
                current_max_id = int(np.max(df['system_idx']))
            #df.drop(columns=['r_inner'], inplace=True)
            df_list.append(df)
            comp_df_list.append(comp_df)
            loop_counts += 1
            if use_pbar:
                pbar.update(gen_stars_chunk)

            neg_missing_stars = np.minimum(missing_stars,0)
            gen_missing_stars = np.maximum(missing_stars,0)

        # combine the results from the different loops
        if len(df_list)==0:
            population_df=pd.DataFrame()
            population_comp_df=pd.DataFrame()
        else:
            population_df = pd.concat(df_list, ignore_index=True)
            if self.mult is not None:
                population_comp_df = pd.concat(comp_df_list, ignore_index=True)
        if self.mult is None:
            population_comp_df = None
        
        # Remove any excess stars
        if (self.lost_mass_option==3) and (len(population_df)>0) and \
            (self.generator.generator_name != 'SpiseaGenerator'):
            population_df = self.remove_stars(population_df, population_comp_df,
                                    neg_missing_stars)
        if len(population_df)>0:
            population_df.loc[:, 'pop'] = self.popid

        #pdb.set_trace()
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
            gg = population_df.groupby(pd.cut(population_df.Dist, radii), observed=False)
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

        logger.log(25, '# Done')
        logger.flush()

        return population_df, population_comp_df

    def generate_stars(self, missing_stars, mass_limit, props, radii=None):
        position = self.population_density.draw_random_positions(missing_stars)
        min_mass = mass_limit

        u, v, w, vr_hc, mu_l, mu_b, vr_lsr = self.do_kinematics(
            position[3], position[4], position[5],
            position[0], position[1], position[2]
            )
        proper_motions = np.column_stack([vr_hc, mu_l, mu_b])
        velocities = np.column_stack([u, v, w, ])

        # generate star at the positions
        star_systems, companions = self.generator.generate_star_at_location(
            position[0:4], props, min_mass, self.max_mass, radii=radii, 
            avg_mass_per_star=self.population_density.average_mass)
        star_systems.loc[:,'x'] = position[0]
        star_systems.loc[:,'y'] = position[1]
        star_systems.loc[:,'z'] = position[2]
        star_systems.loc[:,'Dist'] = position[3]
        star_systems.loc[:,'l'] = position[4]
        star_systems.loc[:,'b'] = position[5]
        star_systems.loc[:,'vr_bc'] = proper_motions[:,0]
        star_systems.loc[:,'mul'] = proper_motions[:,1]
        star_systems.loc[:,'mub'] = proper_motions[:,2]
        star_systems.loc[:,'U'] = velocities[:,0]
        star_systems.loc[:,'V'] = velocities[:,1]
        star_systems.loc[:,'W'] = velocities[:,2]
        star_systems.loc[:,'VR_LSR'] = vr_lsr
        
        if self.obsmag:
            dist_modulus = 5*np.log10(position[3] * 100)
            for band in self.bands:
                star_systems.loc[:,band] += dist_modulus
                if companions is not None:
                    sys_idxs = star_systems['system_idx'].to_numpy()
                    dist_modulus_series = pd.Series(dist_modulus, index=sys_idxs)
                    companions.loc[:,band] += dist_modulus_series[
                            companions['system_idx'].to_numpy()].to_numpy()

        return star_systems, companions

    @staticmethod
    def remove_stars(
            df: pd.DataFrame,
            comp_df: pd.DataFrame,
            missing_stars: np.ndarray,
            ) -> pd.DataFrame:
        """
        Removes stars form data frame if missing_stars is < 0
        """
        if missing_stars<0:
            df = df.sample(n=len(df)+missing_stars, replace=False)
            return df, comp_df[np.isin(comp_df['system_idx'], df['system_idx'])]
        else:
            return df, comp_df

    def check_field(
            self,
            average_imass_per_star: float,
            m_initial: np.ndarray,
            m_evolved: np.ndarray,
            total_mass: float,
            fract_mass_limit: Union[Tuple[float, float], Tuple[np.ndarray, np.ndarray]],
            ) -> np.ndarray:
        """
        Check if stars are missing in a given slice
        Estimates the number of stars that are missing
        Only used for lost_mass_option 3

        Parameters
        ----------
        m_initial : ndarray
        m_evolved : ndarray
        total_mass: float
        fract_mass_limit : ndarray or float
        Returns
        -------
        missing_stars : ndarray
            number of missing stars
        """
        # estimate current initial mass
        #m_in = np.array([np.sum(m_initial[radii_star == r]) for r in radii[:-1]])
        # estimate current evolved mass
        #m_evo = np.array([np.sum(m_evolved[radii_star == r]) for r in radii[:-1]])

        if self.population_density.density_unit in ['number', 'init_mass']:
            return 0

        # option 3.
        # if the sample is large enough otherwise determine the average mass before.
        # otherwise use the determined stars to estimate the average evolved mass:
        # only pick a second time
        average_emass_per_star_mass = (
            fract_mass_limit[0] * fract_mass_limit[1]  # not generated
            + (1 - fract_mass_limit[1]) * np.mean(m_evolved)  # generated
            )
        n_star_expected = total_mass / average_emass_per_star_mass
        total_stars = np.random.poisson(
            n_star_expected * (1 - fract_mass_limit[1]) / self.scale_factor)

        # reduce number of stars by the scale factor
        missing_stars = total_stars - len(m_evolved)
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
        vr_lsr = self.coord_trans.vr_bc_to_vr_lsr(star_l_deg, star_b_deg, vr)

        return u, v, w, vr, mu_l, mu_b, vr_lsr

    def apply_extinction(self, df, comp_df):
        galactic_coordinates = df[['Dist', 'l','b']].to_numpy()

        extinction_in_map, extinction_dict = self.extinction.get_extinctions(
            galactic_coordinates[:, 1],
            galactic_coordinates[:, 2],
            galactic_coordinates[:, 0])

        if self.glbl_params.obsmag:
            for i, band in enumerate(self.bands):
                ext_band = extinction_dict.get(band, 0)
                df.loc[:,band] += ext_band
                if comp_df is not None and (len(comp_df)>0):
                    sys_idxs = df['system_idx'].to_numpy()
                    ext_band_series = pd.Series(ext_band, index=sys_idxs)
                    comp_df.loc[:,band] += ext_band_series[comp_df['system_idx'].to_numpy()].to_numpy()

        df.loc[:,self.extinction.A_or_E_type] = extinction_in_map

        return df, comp_df


