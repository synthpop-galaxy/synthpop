"""
This file contains the Population Class.
Its  task is to handel the generation process for each population
according to the dedicated moduls for
1) Population Density,
2) Initial mass function
3) Age distribution
4) Metallicity distribution
5+6) Isochrone systems and Interpolator
7+8) Extinction Map and Law
9) Kinematic
"""

__all__ = ["Population"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-07-06"
__license__ = "GPLv3"
__version__ = "0.1.0"


# Standard Imports
from posixpath import supports_unicode_filenames
import time
from typing import Tuple, Dict, Set
# Non-Standard Imports
import numpy as np
import pandas

# Local Imports
from . import synthpop_utils as sp_utils
from .position import Position
from .synthpop_utils.synthpop_logging import logger
from .synthpop_utils import coordinates_transformation as coord_trans
from .synthpop_utils import Parameters

# level depended imports based on imported or executed
try:  # check synthpop is imported
    from .. import constants as const
except (ImportError, ValueError):  # if synthpop is executed
    import constants as const
    from modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from modules.evolution import EvolutionIsochrones, EvolutionInterpolator, \
        CombineEvolution, MUST_HAVE_COLUMNS
    from modules.age import Age
    from modules.initial_mass_function import InitialMassFunction
    from modules.kinematics import Kinematics
    from modules.metallicity import Metallicity
    from modules.population_density import PopulationDensity

else:  # continue import when if synthpop is imported
    from ..modules.extinction import ExtinctionLaw, ExtinctionMap, CombineExtinction
    from ..modules.evolution import EvolutionIsochrones, EvolutionInterpolator, \
        CombineEvolution, MUST_HAVE_COLUMNS
    from ..modules.age import Age
    from ..modules.initial_mass_function import InitialMassFunction
    from ..modules.kinematics import Kinematics
    from ..modules.metallicity import Metallicity
    from ..modules.population_density import PopulationDensity


class Population:
    """
    Population subclass for model class.
    Requires a population parameter file to be initialized.

    ...

    Attributes
    ----------
    name : string
    l_deg : float
    b_deg : float
    solid_angle : float
    coordinate_flag : bool
    popid : int
    params : dict

    age : age.Age class
    imf : initial_mass_function.InitialMassFunction subclass
    kinematics : kinematics.Kinetamics subclass
    metallicity : metallicity.Metallicity subclass
    population_density : population_density.PopulationDensity subclass
    position : position.Position class

    Methods
    -------
    __init__(self,population_file,population_index,isochrones,bands,**kwargs) : None
        initialize the Population class
    assign_subclasses(self, subclasses) : None
        convenience function to assign subfuctions to Population upon initialization
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
    generate_star_position(self,r_inner,r_outer, N) : float x 8
        generate the position, velocity, distance, coordinates for N stars
    evolve_field(self,data) : float x 4
        estimates the evolved parameter for each star in the field
    """

    # still deciding on what goes in init...
    # we could also make an update function to reset values like coordinates and solid angle
    # new args ideas : population_file, field class, evolution class, extinction class
    def __init__(
            self, population_file: str, population_index: int, glbl_params: Parameters,
            **positional_kwargs
            ):
        # machinery to import parameter file for a population only given filename
        self.glbl_params = glbl_params

        # reduces the number of stars by a given scale_factor
        # i.e. if set to 100 only 1% of the stars are generated
        self.scale_factor = getattr(glbl_params, "scale_factor", 1)

        self.pop_params = sp_utils.json_loader(population_file)
        self.name = self.pop_params['name']
        self.popid = population_index

        logger.create_info_subsection(f"Population {self.popid};  {self.name}")
        # independent index issued by model class

        logger.info(f"Initialize Population {self.popid} ({self.name}) from {population_file}")
        logger.debug("pop_params = {")
        for key,val  in self.pop_params.items():
            if isinstance(val, str): val=f"'{val}'"
            logger.debug(f'    "{key}": {val},')
        logger.debug("    }")

            # default pointing is Galactic anti-center
        self.l_deg = 180
        self.b_deg = 0
        self.solid_angle_sr = 1e-12
        self.position = Position(self.l_deg, self.b_deg, self.solid_angle_sr,
            **{**self.glbl_params.window_type, **positional_kwargs}
        )
        self.coordinate_flag = False

        # read min and maximum mass from parameters
        self.min_mass = self.glbl_params.mass_lims['min_mass']
        self.max_mass = self.glbl_params.mass_lims['max_mass']
        self.skip_lowmass_stars = self.glbl_params.skip_lowmass_stars
        # get maximum distance and step_size
        self.max_distance = self.pop_params.get('max_distance', self.glbl_params.max_distance)
        self.step_size = self.pop_params.get('distance_step_size',
            self.glbl_params.distance_step_size)

        self.N_mc_totmass = self.glbl_params.N_mc_totmass

        # check to see if the filename matches the population name in file,
        # just a check for errors if copy/paste-ing
        # if self.pop_params['name'].lower() != population_file.split('.')[0].lower():
        if population_file.split('.')[0].lower() != population_file.split('.')[0].lower():
            logger.warning(f"Input population filename ({population_file.split('.')[0].lower()})"
                           " does not match name in file ({self.pop_params['name'].lower()}). ")

        # star number correction
        self.lost_mass_option = self.pop_params.get("lost_mass_option",
            self.glbl_params.lost_mass_option)

        self.N_av_mass = self.pop_params.get("N_av_mass", self.glbl_params.N_av_mass)

        # placeholder for profiles etc. initialized via the assign_subclasses
        self.evolution = None
        self.extinction = None
        self.age = None
        self.imf = None
        self.kinematics = None
        self.metallicity = None
        self.assign_subclasses()
        if isinstance(self.evolution, list):
            self.bands = list(self.evolution[0].bands)
        else:
            self.bands = list(self.evolution.bands)
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

    def assign_subclasses(self):
        """initialization of all the subclass"""

        # get evolution class
        evolution_class = None
        if hasattr(self.glbl_params, "evolution_class"):
            evolution_class = self.glbl_params.evolution_class
        evolution_class = getattr(self.glbl_params, "evolution_class")
        if evolution_class is None:
            # log N
            if "evolution_kwargs" in self.pop_params:
                logger.debug("read evolution class from Population file " )
                evolution_class= self.pop_params.get("evolution_kwargs")
            else:
                msg = "evolution_kwargs is neither defined " \
                      "via the config file nor via the population"
                logger.critical(msg)
                raise ModuleNotFoundError(msg)
        else:
            logger.debug("read evolution class from config file ")

        if isinstance(evolution_class, dict):
            evolution_class = [evolution_class]

        evolutions = []
        # go through  all mass ranges
        # this is to combine different isochrone systems or interpolator
        ev_init = []
        for i, ev_kwargs in enumerate(evolution_class):
            # get the Isochrone system
            Isochrone_System = sp_utils.get_subclass(
                EvolutionIsochrones, ev_kwargs,
                initialize=False, population_file=self.pop_params['__filename'])

            # check if an interpolator is specified
            if 'interpolator' in ev_kwargs.keys():
                Interpolator = sp_utils.get_subclass(
                    EvolutionInterpolator, {"name": ev_kwargs['interpolator']},
                    initialize=False, population_file=self.pop_params['__filename'])
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
                # convert dict into list of tuples
                columns.extend(self.glbl_params.chosen_bands.items())
            else:
                columns.extend(self.glbl_params.chosen_bands)
            ev_kwargs.update({"columns": columns})

            ev_init.append(ev_kwargs)
            evolution_i = Evo(iso_kwargs=ev_kwargs, int_kwargs=ev_kwargs)
            val = ev_kwargs.get("min_mass", evolution_i.min_mass)
            if val < evolution_i.min_mass:
                msg = f'min mass ({val:.3f}) in the evolution keywords ' \
                      f'is lower than the masses in the Isochrone grid ({evolution_i.min_mass:.3f})'
                logger.critical(msg)
                raise ValueError(msg)

            evolution_i.min_mass = val

            evolution_i.max_mass = ev_kwargs.get("max_mass", evolution_i.max_mass)
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
            self.evolution = evolutions[0]
        else:
            self.evolution = evolutions

        # get Extinctions
        ExtMap = sp_utils.get_subclass(ExtinctionMap, self.glbl_params.extinction_map_kwargs,
            initialize=False, population_file=self.pop_params['__filename'])

        logger.debug(
            [f'initialize Extinction Map with keywords {self.glbl_params.extinction_map_kwargs}',
                ''])
        logger.log(15, f'"extinction_map" :  {self.glbl_params.extinction_map_kwargs}')

        # check if multiple extinction laws are provided
        if isinstance(self.glbl_params.extinction_law_kwargs, dict):
            ExtLaw = sp_utils.get_subclass(ExtinctionLaw, self.glbl_params.extinction_law_kwargs,
                 initialize=False, population_file=self.pop_params['__filename'])
        else:
            # collect all extinction laws in a common list
            ExtLaw = [
                sp_utils.get_subclass(ExtinctionLaw, ext_law, initialize=False,
                    population_file=self.pop_params['__filename'])
                for ext_law in self.glbl_params.extinction_law_kwargs]

        # combine Extinction and Extinction laws in a combined Class
        Extinction = CombineExtinction(ext_map=ExtMap, ext_law=ExtLaw,
            which_extinction_law=self.glbl_params.extinction_law_kwargs)
        # initialize extinction

        logger.log(15, '"extinction_law_kwargs" : [')
        for ext_law_kwarg in self.glbl_params.extinction_law_kwargs:
            msg = f'    {ext_law_kwarg},'.replace("'", '"')
            msg = msg.replace('False', 'false').replace('True', 'true')
            msg = msg.replace('None', 'null')
            logger.debug(f'initialize Extinction law with keywords {ext_law_kwarg}')
            logger.log(15, msg)
        logger.log(15, '   ]')

        self.extinction = Extinction(ext_map_kwargs=self.glbl_params.extinction_map_kwargs,
            ext_law_kwargs=self.glbl_params.extinction_law_kwargs)
        self.extinction.set_R_V(self.glbl_params.R_V)
        # assign age
        logger.debug("%s : using Age subclass '%s'", self.name,
            self.pop_params['age_func_kwargs']['name'])
        self.age = sp_utils.get_subclass(
            Age, self.pop_params['age_func_kwargs'],
            keyword_name='age_func_kwargs',
            population_file=self.pop_params['__filename'])

        logger.debug("%s : using InitialMassFunction subclass '%s'", self.name,
            self.pop_params['imf_func_kwargs']['name'])
        self.imf = sp_utils.get_subclass(
            InitialMassFunction, self.pop_params['imf_func_kwargs'],
            keyword_name='imf_func_kwargs',
            population_file=self.pop_params['__filename'])

        logger.debug("%s : using Kinematics subclass '%s'", self.name,
            self.pop_params['kinematics_func_kwargs']['name'])
        self.kinematics = sp_utils.get_subclass(
            Kinematics, self.pop_params['kinematics_func_kwargs'],
            keyword_name='metallicity_func_kwargs',
            population_file=self.pop_params['__filename'])

        logger.debug("%s : using Metallicity subclass '%s'", self.name,
            self.pop_params['metallicity_func_kwargs']['name'])
        self.metallicity = sp_utils.get_subclass(
            Metallicity, self.pop_params['metallicity_func_kwargs'],
            keyword_name='metallicity_func_kwargs',
            population_file=self.pop_params['__filename'])

        logger.debug("%s : using PopulationDensity subclass '%s'", self.name,
            self.pop_params['population_density_kwargs']['name'])
        self.population_density = sp_utils.get_subclass(
            PopulationDensity, self.pop_params['population_density_kwargs'],
            keyword_name='population_density_kwargs',
            population_file=self.pop_params['__filename'])

    def set_position(self, l_deg: float, b_deg: float, solid_angle: float, solid_angle_unit: str,
            **position_kwargs):
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
            self.solid_angle_sr = solid_angle * (np.pi/180)**2
        else:
            raise ValueError(f"solid angle unit '{solid_angle_unit}' is not suported")

        self.coordinate_flag = True
        self.position.update(l_deg, b_deg, self.solid_angle_sr,
            **{**self.glbl_params.window_type,**position_kwargs})
        self.extinction.update_line_of_sight(l_deg, b_deg)
        # flag to mark that coordinates have been set
        logger.debug(
            f"{self.name} : position is set to "
            f"{self.l_deg: .3f}, {self.b_deg: .3f} with solid angle"
            f"{self.solid_angle_sr:.3e} sr ({solid_angle:.3e} {solid_angle_unit})")

    def mc_totmass(self, r_inner: float, r_outer: float, n_picks: int = 1000) -> float:
        """
        Monte Carlo integration of the total mass in a slice,
        given near and far bounding radii r_inner and r_outer
        Want to add some catch for large error in the mass,
        then increase value of N
        Recall that for a Monte Carlo integration
        Q_N = Volume * (1/N) * sum_N f(x) = V*<f>
        the expected error in Q_N decreases as 1/sqrt(N)

        Parameters
        ----------
        r_inner : float
            inner radius of the slice
        r_outer : float
            outer radius of the slice
        n_picks : int
            number of picks foor the  integration

        Returns
        -------
        total_mass : float
            total mass in the slice
        """
        # calculate volume of slice, needs to be kpc-3, r_inner and r_outer in kpc
        volume = (1 / 3) * self.solid_angle_sr * (r_outer ** 3 - r_inner ** 3)

        # MC draws of density
        d, _, _, _, lstar_deg, bstar_deg = self.position.draw_random_point_in_slice(r_inner,
            r_outer, n_picks)
        r, phi_rad, z = coord_trans.dlb_to_rphiz(d, lstar_deg, bstar_deg)
        mean_density = np.mean(self.population_density.density(r, phi_rad, z))
        # then find mass from volume*<density>
        total_mass = volume * mean_density

        return total_mass

    def central_totmass(self, r_inner: float, r_outer: float) -> float:
        """
        Monte Carlo integration of the total mass in a slice,
        given near and far bounding radii r_inner and r_outer
        Want to add some catch for large error in the mass, then increase value of N
        Recall that for a Monte Carlo integration
        Q_N = Volume * (1/N) * sum_N f(x) = V*<f>
        the expected error in Q_N decreases as 1/sqrt(N)
        """

        # calculate volume of slice, needs to be kpc-3, r_inner and r_outer in kpc
        volume = (1 / 3) * self.solid_angle_sr * ((r_outer) ** 3 - (r_inner) ** 3)
        r, phi_rad, z = coordinates.dlb_to_rphiz((r_outer + r_inner) / 2, self.l_deg, self.b_deg)
        density = self.population_density.density(r, phi_rad, z)

        totmass = volume * density

        return totmass

    def estimate_field(self, **kwargs) -> Tuple[float, float]:
        """
        estimate the field in the current pointing
        """

        if not self.coordinate_flag:
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
            av_mass_corr = self.estimate_evolved_average_mass(n_stars=self.N_av_mass,
                use_save_value=use_save_av_mass)
            logger.debug('%s : average mass loss: %.1f%s', self.name, av_mass_corr * 100, '%')
        else:
            av_mass_corr = self.pop_params.get("av_mass_corr",
                1 / self.pop_params.get('n_star_corr', 1.4))
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

    def estimate_evolved_average_mass(
            self,
            n_stars: int = 10000,
            use_save_value: bool = False,
            **kwargs
            ) -> float:
        """
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
            return self.av_mass_corr
        data = self.generate_stars_in_slice(0, self.max_distance, n_stars)
        m_initial, age, metallicity, _, _ = data[:, [-3, -2, -1, 0, 1]].T

        _, s_props, _, _ , not_evolved = self.multi_get_evolved_props(
            m_initial, metallicity, age, {const.REQ_ISO_PROPS[0], self.glbl_params.maglim[0]})

        mass_evolved = s_props[const.REQ_ISO_PROPS[0]]
        mass_evolved[not_evolved] = m_initial[not_evolved]
        average_m_initial = self.imf.average_mass(min_mass=self.min_mass, max_mass=self.max_mass)

        av_mass_corr = np.mean(mass_evolved) / average_m_initial

        self.av_mass_corr = av_mass_corr

        return av_mass_corr

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
        if not self.coordinate_flag:
            msg = ("coordinates where not set."
                   "You must run set_position(l_deg, b_deg, solid_angle_sr)"
                   "before running 'estimate_field' or 'generate_field'")
            logger.critical(msg)
            raise AttributeError(msg)

        # array of radii to use
        logger.create_info_subsection(f"Population {self.popid};  {self.name}")
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

        # estimate the mass/numbers in each slice
        mass_per_slice = np.array([self.mc_totmass(radii_inner, radii_outer, self.N_mc_totmass)
            for radii_inner, radii_outer in zip(radii, radii[1:])])

        # get average_mass from imf
        average_imass_per_star = self.imf.average_mass(min_mass=self.min_mass,
            max_mass=self.max_mass)

        ################################################################
        #          Estimate mass loss correction factor                #
        ################################################################

        if self.population_density.density_unit == 'init_mass':
            av_mass_corr = 1
        # check which correction should be performed

        elif self.lost_mass_option == 1:
            # estimate the av_mass_corr in the estimate field
            # and use it for different positions
            av_mass_corr = self.estimate_evolved_average_mass(
                n_stars=self.N_av_mass, use_save_value=True
                )

        elif self.lost_mass_option == 2:
            # estimate the av_mass_corr for each field
            av_mass_corr = self.estimate_evolved_average_mass(
                n_stars=self.N_av_mass, use_save_value=True
                )

        elif self.lost_mass_option == 3:
            # estimate the av_mass_corr for each field,
            # but uses the data if the number of stars N_av_mass == True
            # will be estimated in the check_field
            # otherwise
            n_star_expected = mass_per_slice / average_imass_per_star
            av_mass_corr = 1
            if np.sum(n_star_expected) < self.N_av_mass:
                n_star_expected *= self.N_av_mass / np.sum(n_star_expected)

        elif self.lost_mass_option == 6:
            # use the n_star_corr keyword specified in the population_config_file
            av_mass_corr = self.pop_params.get("av_mass_corr",
                1 / self.pop_params.get('n_star_corr', 1.4))
        else:
            av_mass_corr = 1

        ti = time.time()  # start timer

        ################################################################
        #   Translate density into number of generated stars           #
        ################################################################
        if self.population_density.density_unit == "number":
            logger.info("# From density profile (number density)")
            n_star_expected = mass_per_slice  # density returns n_star_expected
            mass_per_slice = n_star_expected * average_imass_per_star * av_mass_corr

        else:
            logger.info("# From density profile")
            n_star_expected = mass_per_slice / (average_imass_per_star * av_mass_corr)

        distribution["distance_distribution"] = np.vstack(
            [(radii[1:] + radii[:-1]) / 2, n_star_expected]).T
        distribution["distance_distribution_comment"] = \
            "pairs of distances in [kpc] and number of stars in the slice"

        mass_limit = np.ones(radii[:-1].shape) * self.min_mass  # new min mass for each
        frac_lowmass = (0, 0)  # average mass/ fraction of stars

        if (self.glbl_params.maglim[1] > 50) \
                or (not self.skip_lowmass_stars) \
                or (not self.glbl_params.obsmag):
            pass
        elif self.skip_lowmass_stars and isinstance(self.evolution, list):
            logger.warning("skip_lowmass_stars is not implemented for multiple isochrones")

        elif self.skip_lowmass_stars:
            logger.info(f"{self.name} : estimate minimum mass for magnitude limit")
            max_age = self.age.get_maximum_age()
            mass_limit = self.evolution.get_mass_min(self.glbl_params.maglim[0],
                self.glbl_params.maglim[1], radii[:-1], max_age)
            mass_limit = np.maximum(mass_limit, self.min_mass)
            mass_limit = np.minimum(mass_limit, self.max_mass)
            frac_lowmass = (
                np.array([self.imf.average_mass(self.min_mass, mm) for mm in mass_limit]),
                (self.imf.F_imf(mass_limit) - self.imf.F_imf(self.min_mass)) /
                (self.imf.F_imf(self.max_mass) - self.imf.F_imf(self.min_mass))
                )  # average_mass, fraction of stars

        # reduce number of stars by the scale factor and fract_above_min_mass
        total_stars = np.random.poisson(n_star_expected * (1 - frac_lowmass[1]) / self.scale_factor)

        expected_total_imass = n_star_expected * average_imass_per_star

        logger.info(f"expected_total_Imass = {np.sum(expected_total_imass):.4f}")
        logger.info(f"expected_total_Emass = {np.sum(mass_per_slice):.4f}")
        logger.info(f"average_Imass_per_star = {average_imass_per_star:.4f}")
        logger.info(f"mass_loss_correction = {np.sum(av_mass_corr):.4f}")
        logger.info(f"n_expected_stars = {np.sum(n_star_expected):.4f}")
        if self.skip_lowmass_stars:
            logger.info(f"without_lm_stars = {np.sum(n_star_expected * (1 - frac_lowmass[1])):.4f}")
        logger.debug(f"{self.name} : Lost mass option: %s", self.lost_mass_option)
        logger.debug(f"{self.name} : Generate ~{sum(total_stars)} stars")
        if not self.glbl_params.kinematics_at_the_end:
            logger.info("# Determine velocities when position are generated ")

        ################################################################
        #                  Generate Stars                              #
        ################################################################
        missing_stars = total_stars
        loop_counts = 0
        logger.debug("generate stellar properties")
        while any(missing_stars > 0):
            pop_list = []
            for radii_inner, radii_outer, n_stars, min_mass \
                    in zip(radii, radii[1:], missing_stars, mass_limit):
                # Generate stars in slice
                if n_stars == 0:
                    continue
                # initial stellar parameters

                init_props = self.generate_stars_in_slice(radii_inner, radii_outer, n_stars,
                    min_mass=min_mass, max_mass=self.max_mass)
                pop_list.append(init_props)
            # combine all slices in a common array of data
            pop_array = np.vstack(pop_list)

            # Evolve the field at once
            props, extinction_in_map, mags, user_props, final_phase_flag = \
                self.evolve_field(pop_array)
            m_evolved, props = props[0], props[1:]

            # check field, eg. does the density matches the expectations
            missing_stars = self.check_field(radii, average_imass_per_star,
                pop_array[:, -3], m_evolved, pop_array[:, 0], loop_counts,
                mass_per_slice, frac_lowmass)
            # Convert Table to pd.DataFrame
            # .T is used to split the pop_array along the columns and not along rows
            # "pop", "iMass", "age",  "Fe/H_initial", "Mass","In_Final_Phase", "Dist", "l", "b",
            #         "mul", "mub", "vr",  "x", "y", "z",  "U", "V", "W", "VR_LSR", "ExtinctionInMap"


            df = pandas.DataFrame(np.vstack(
                [self.popid * np.ones(len(pop_array)),  # pop,
                *pop_array[:, [14, 15, 16]].T, m_evolved,  final_phase_flag, # iMass, age, Fe/H,  Mass, In_Final_Phase
                *pop_array[:, [1, 5, 6, 10, 11, 12]].T,  # Dist, l, b, Vr, mu_l, mu_b,
                *pop_array[:, [2, 3, 4, 7, 8, 9, 13]].T,  #  x, y, z,  U, V, W, VR_LSR
                extinction_in_map,  # extinction_in_map
                *props, *user_props, *mags  # all specified props and magnitudes
                ]
                ).T, columns=headers)

            # remove stars if we have overestimated the number of stars
            # might happen for mass loss option 3 and 4:
            if any(missing_stars < 0):
                df = self.remove_stars(df, pop_array[:, 0], missing_stars, radii)
                missing_stars = np.maximum(missing_stars, 0)
            # add to previous drawn data
            df_list.append(df)
            loop_counts += 1
        if len(df_list) == 0:
            population_df = pandas.DataFrame(columns=headers, dtype=float)
        else:
            population_df = pandas.concat(df_list, ignore_index=True)
        to = time.time()  # end timer

        logger.debug("%s : Generated %i stars in %i cycles (%.2fs)", self.name, len(population_df),
            loop_counts, to - ti)
        logger.info('# From Generated Field:')

        logger.info(f'generated_stars = {len(population_df)}')
        if len(population_df) != 0:

            logger.info(f'generated_total_iMass = {population_df["iMass"].sum():.4f}')
            if self.skip_lowmass_stars:
                gg = population_df.groupby(pandas.cut(population_df.Dist, radii))

                im_incl = (gg["iMass"].sum()
                           + gg.size() * frac_lowmass[0] * frac_lowmass[1] / (
                                   1 - frac_lowmass[1])).sum()

                logger.info(f'generated_total_iMass_incl_lowmass = {im_incl.sum():.4f}')

            logger.info(f'generated_total_eMass = {population_df["Mass"].sum():.4f}')
            if self.skip_lowmass_stars:
                em_incl = (gg["Mass"].sum()
                           + gg.size() * frac_lowmass[0] * frac_lowmass[1] / (
                                   1 - frac_lowmass[1])).sum()
                logger.info(f'generated_total_eMass_incl_lowmass = {em_incl:.4f}')

            logger.info(f'det_mass_loss_corr = '
                        f'{population_df["Mass"].sum() / population_df["iMass"].sum():.4f}')
            if self.skip_lowmass_stars:
                logger.info(f'det_mass_loss_corr_incl_lowmass = {em_incl / im_incl}')

            logger.debug(f'average_mass_per_star = {population_df["Mass"].mean():.4f}')
            if self.skip_lowmass_stars:
                mean_mass = em_incl * ((1 - frac_lowmass[1]) / gg.size()).sum()
                logger.debug(f'average_mass_per_star_incl_lowmass = {mean_mass:.4f}')

        if (self.glbl_params.maglim[-1] != 'keep'):
            criteria = population_df[self.glbl_params.maglim[0]] < self.glbl_params.maglim[1]
        else:
            criteria = None

        sp_utils.log_basic_statistics(population_df, f"stats_{self.name}", criteria)
        logger.log(25,'# Done')
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
            loop_counts: int,
            mass_per_slice: np.ndarray,
            fract_mass_limit: np.ndarray or float = 1.,
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
            n_star_expected = mass_per_slice / (average_emass_per_star_mass)
            total_stars = np.random.poisson(
                n_star_expected * (1 - fract_mass_limit[1]) / self.scale_factor)

            # reduce number of stars by the scale factor
            exist = np.array([np.sum(radii_star == r) for r in radii[:-1]])
            missing_stars = total_stars - exist
            return missing_stars

        else:
            return np.zeros(1)

    def generate_stars_in_slice(
            self, r_inner: float, r_outer: float, n_stars: None or int = None,
            min_mass: None or float = None, max_mass: None or float = None
            ):
        """
        Generates a the initial parameters of N stars between r_inner and r_outer

        Parameters
        ----------
        r_inner : float:
            inner radius of slice
        r_outer : float
            outer radius of slice
        n_stars : int, None, optional
            Number of stars generated
            if None, a single  floats will be return for each output

        Returns
        -------
        r_in : float,nd_array
            r_inner to identify the slice the stars belong to
        dist : float, ndarray [kpc]
            distance
        x, y, z, : float, ndarray [kpc]
            galactocentric cartesian coordinates
        star_l_deg, star_b_deg : float, ndarray [degree]
            galactic longitude and latitude
        u, v, w : float, ndarray [km/s]
            velocitys
        vr : float, ndarray [km/s]
            radial velocity
        mu_l, mu_b : float, ndarray [mas/yr]
            proper motion in galactic coordinates
        m_initial : float, ndarray [Msun]
            inital mass in solar masses
        age : float, ndarray [Gyr]
            age of the star
        met :  float, ndarray
            initial metalicitiy of the stars
        """
        # pos_data = dist, x, y, z, star_l_deg, star_b_deg, u, v, w, vr, mu_l, mu_b
        pos_data = self.generate_star_position(r_inner, r_outer, n_stars)

        # draw initial properties from distribution
        if min_mass is None:
            min_mass = self.min_mass
        if max_mass is None:
            max_mass = self.max_mass
        m_initial = self.imf.draw_random_mass(
            N=n_stars, min_mass=min_mass, max_mass=max_mass)
        age = self.age.draw_random_age(N=n_stars)
        metallicity = self.metallicity.draw_random_metallicity(
            N=n_stars, x=pos_data[1], y=pos_data[2], z=pos_data[3])

        if n_stars is None:
            return (r_inner, *pos_data, m_initial, age, metallicity)
        else:
            return np.vstack([np.ones(n_stars) * r_inner, *pos_data, m_initial, age, metallicity]).T

    def generate_star_position(self, r_inner: float, r_outer: float, n_stars: int = 1) \
            -> Tuple[np.ndarray,]:
        """
        generates the star position and proper motion for N stars within a slice

        Parameters
        ----------
        r_inner : float:
            inner radius of slice
        r_outer : float
            outer radius of slice
        n_stars : int, None, optional
            Number of stars generated
            if None, a single  floats will be return for each output

        Returns
        -------
        dist : float, ndarray [kpc]
            distance
        x, y, z, : float, ndarray [kpc]
            galactocentric cartesian coordinates
        star_l_deg, star_b_deg : float, ndarray [degree]
            galactic longitude and latitude
        u, v, w : float, ndarray [km/s]
            velocitys
        vr : float, ndarray [km/s]
            radial velocity
        mu_l, mu_b : float, ndarray [mas/yr]
            proper motion in galactic coordinates
        """

        if n_stars <= 0:
            return

        # for rejection sampling, we will use a uniform distribution with
        # amplitude N the density of the center of the slice
        # then, we will generate random position, find that locations density,
        # and compare to uniform draw until repeat
        # rejection based on 1.25*(a randomly drawn density)
        _, x_reject, y_reject, z_reject, _, _ = (
            self.position.draw_random_point_in_slice(r_inner, r_outer, n_stars))
        rejection_density = 1.25 * self.population_density.density(
            np.sqrt(x_reject ** 2 + y_reject ** 2), np.arctan2(y_reject, x_reject), z_reject)

        position_not_selected = np.ones(n_stars, bool)
        data = np.zeros((n_stars, 6))
        while any(position_not_selected):
            randDensity = np.random.uniform(0, rejection_density[position_not_selected])
            data_step = np.vstack(
                self.position.draw_random_point_in_slice(r_inner, r_outer,
                    n_stars=np.sum(position_not_selected))
                ).T
            drawDensity = self.population_density.density(
                np.sqrt(data_step[:, 1] ** 2 + data_step[:, 2] ** 2),
                np.arctan2(data_step[:, 2], data_step[:, 1]), data_step[:, 3])
            # check which one can be used
            not_take = randDensity > drawDensity
            data[position_not_selected] = data_step
            position_not_selected[position_not_selected] = not_take

        # unpack data
        dist, x, y, z, star_l_deg, star_b_deg = data.T

        if self.glbl_params.kinematics_at_the_end:
            u, v, w, vr_hc, mu_l, mu_b, vr_lsr = np.zeros((7, n_stars))
        else:
            u, v, w, vr_hc, mu_l, mu_b, vr_lsr = self.do_kinematics(
                dist, star_l_deg, star_b_deg, x, y, z)

        return dist, x, y, z, star_l_deg, star_b_deg, u, v, w, vr_hc, mu_l, mu_b, vr_lsr

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

        # draw random velocities,
        u, v, w = self.kinematics.draw_random_velocity(x, y, z, **kwargs)

        # translate u, v, w into mu_l and mu_b_vr
        vr, mu_l, mu_b = coord_trans.uvw_to_vrmulb(star_l_deg, star_b_deg, dist, u, v, w)

        # correct for motion of the sun
        # following Beaulieu et al. (2000)
        vr_lsr = (vr
              + const.V_LSR * np.sin(star_l_deg * np.pi / 180)
              * np.sin(star_b_deg * np.pi / 180)
              + const.V_PEC * (
                      np.sin(star_b_deg * np.pi / 180) * np.sin(const.B_APX_DEG * np.pi / 180)
                      + np.cos(star_b_deg * np.pi / 180) * np.cos(const.B_APX_DEG * np.pi / 180)
                      * np.cos(star_l_deg * np.pi / 180 - const.L_APX_DEG * np.pi / 180)
              )
              )
        return u, v, w, vr, mu_l, mu_b, vr_lsr

    def evolve_field(self, data: np.ndarray):
        """
        Starts the evolution process and performs the edge handling
        translates from absolute to relative coordinates
        """
        # combine columns which should be interpolated
        # use set to avoid duplicates
        props = set(const.REQ_ISO_PROPS + self.glbl_params.opt_iso_props + self.bands)

        # extract data from
        m_initial, age, metallicity, radii_inner, distance, l_deg, b_deg = \
            data[:, [-3, -2, -1, 0, 1, 5, 6]].T

        # default values
        # default_props[:,0] is the stellar mass
        default_props = np.zeros((len(const.REQ_ISO_PROPS), len(m_initial)))

        mags = 9999 * np.ones((len(self.bands), len(m_initial)))
        extinction_in_map = np.zeros(len(m_initial))
        user_props = np.zeros((len(self.glbl_params.opt_iso_props), len(m_initial)))

        init_mag, s_props, final_phase_flag,  inside_grid, not_evolved = self.multi_get_evolved_props(
            m_initial, metallicity, age, props)
        if self.glbl_params.obsmag:
            init_mag[inside_grid] += 5 * np.log10(distance[inside_grid] * 100)

        # sort interpolated properties
        for i, band in enumerate(self.bands):
            mags[i] = s_props[band]
        for i, key in enumerate(const.REQ_ISO_PROPS):
            default_props[i] = s_props[key]
        default_props[1:, np.logical_not(inside_grid)] = np.nan
        default_props[1:, not_evolved] = np.nan
        default_props[0, not_evolved] = m_initial[not_evolved]

        for i, key in enumerate(self.glbl_params.opt_iso_props):
            user_props[i] = s_props[key]
        user_props[:, np.logical_not(inside_grid)] = np.nan
        user_props[:, not_evolved] = np.nan

        # convert magnitude in observed magnitudes
        for ri in np.unique(radii_inner):
            # select stars in current slice
            current_slice = radii_inner == ri

            # update Extinctions
            self.extinction.update_extinction_in_map(radius=ri)

            # get extinctions for stars in current slice
            extinction_in_map[current_slice], extinction_dict = self.extinction.get_extinctions(
                l_deg[current_slice],
                b_deg[current_slice],
                distance[current_slice])

            # translate magnitudes to observed magnitudes
            if self.glbl_params.obsmag:
                # check if magnitude is less than magnitude limits
                ext_mag = extinction_dict.get(self.glbl_params.maglim[0], 0)
                init_mag[current_slice] += ext_mag

                for i, band in enumerate(self.bands):
                    mags[i, current_slice] += 5 * np.log10(distance[current_slice] * 100) \
                                              + extinction_dict.get(band, 0)

        mag_le_limit = init_mag < self.glbl_params.maglim[1]
        mags[:, np.logical_not(mag_le_limit)] = np.nan

        return default_props, extinction_in_map, mags, user_props, final_phase_flag

    def multi_get_evolved_props(
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
        logger.debug("Start evolving field")
        ti=time.time()
        s_track = {p: np.ones(len(m_init)) * np.nan for p in props}
        mag = np.nan * np.ones(len(m_init))
        inside_grid = np.ones(len(m_init), bool)
        in_final_phase = np.ones(len(m_init), bool)
        not_performed = np.ones(len(m_init), bool)
        evolution_list = self.evolution if isinstance(self.evolution, list) else (self.evolution,)
        for i, evolution_i in enumerate(evolution_list):
            # check if evolution has an atribute which says if numpy arrays can be used
            if hasattr(evolution_i, 'accept_np_arrays'):
                accept_np_arrays = evolution_i.accept_np_arrays
            else:
                accept_np_arrays = True

            # find the stars which falles into the mass range of the current evolution class
            if i != len(evolution_list) - 1:
                which = np.where(not_performed & (m_init > evolution_i.min_mass) & (
                        m_init < evolution_i.max_mass))[0]
            else:
                which = np.where(not_performed & (m_init > evolution_i.min_mass))[0]

            # check if there are any stars for this step
            if len(which) == 0:
                continue
            # loop over bunches of at most 250k to reduce memory usage
            if len(which) > self.glbl_params.chunk_size * 6 / 5:
                chunk_size = self.glbl_params.chunk_size
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
                    s_props_i, inside_grid_i , in_final_phase_i = evolution_i.get_evolved_props(
                        m_init[which2], met[which2], age[which2], props, **kwargs)

                else:
                    # This can be used if self.evolution.get_evolved_props
                    # can not handle multiple stars and numpy array:
                    m_initial = m_init[which2]
                    metallicity = met[which2]
                    age2 = age[which2]
                    s_props_array, inside_grid_i , in_final_phase_i= np.array([
                        self.evolution.get_evolved_props(
                            m_initial[i], metallicity[i], age2[i] * 1e9, props, **kwargs)
                        for i in range(len(m_initial))
                        ]).T
                    # s_props needs to be transformed from an array of dictionaries
                    # into a dictionary of numpy arrays
                    s_props_i = {
                        key: np.array([i[key] for i in s_props_array])
                        for key in s_props_array[0].keys()}

                # update results to data array
                mag_i = s_props_i.get(self.glbl_params.maglim[0])
                mag[which2] = mag_i
                inside_grid[which2] = inside_grid_i
                in_final_phase[which2] = in_final_phase_i
                for key in s_track.keys():
                    s_track[key][which2] = s_props_i[key]
                # update the list of not performed stars

                count_c += len(which2)
                if use_chunks:
                    print("\r", count_c, "/", len(which), end='')

            # update the list of not performed stars
            not_performed[which] = False
            print('')
            # check if anything left to do
            if not any(not_performed):
                break
        logger.debug(f"used time = {time.time()-ti:.2f}s")

        return mag, s_track, in_final_phase, inside_grid, not_performed