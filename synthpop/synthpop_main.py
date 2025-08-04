"""
SynthPop is a modular framework to generate synthetic galaxy population models.
For usage see README.md!

This file contains the main SynthPop class and main function.
Which handles the setting of synthpop, data collection
from the different populations and saving process.
The generation process for each population is performed by
the Population class defined in population.py.
"""

__all__ = ["SynthPop", "main"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__data__ = "2023-01-09"
__license__ = "GPLv3"
__version__ = "1.0.0"

# Standard Imports
import os
import time
import datetime
import glob
import importlib
from typing import Tuple, Dict, Iterator
# Non-Standard Imports
import pandas
import numpy as np

# check if astropy is installed
if importlib.util.find_spec("astropy") is not None:
    import astropy.table as astrotable

else:
    pass
    # astropy is only needed if the output format is either VoTable or FITS-table

# Local Imports
try:
    from . import constants as const

except (ImportError, ValueError) as e:
    import constants as const
    import synthpop_utils as sp_utils
    from modules.post_processing import PostProcessing
    from population import Population
    from synthpop_utils.synthpop_logging import logger

else:
    from .modules.post_processing import PostProcessing
    from . import synthpop_utils as sp_utils
    from .population import Population
    from .synthpop_utils.synthpop_logging import logger


class SynthPop:
    """
    Model class for generating catalogs of stars.

    Attributes
    ----------
    model_name : str
        name of the model directory
    population_files : list<str>
        list of population.json files detected in the model directory
    population : list<Population>
        list of the initialized population objects
    solid_angle : float [sr]
        size of the cone
    min_mass, max_mass : float [Msun]
        minimum and maximum mass of the Star generation
    filename_base : str
        file name to store the results including the  path and  without an  extension

    Methods
    -------
    __init__() : None
        initialize the SynthPop class
    init_populations() : None
        collects the population configuration and initializes the populations
    update_location(l_deg: float, b_deg: float, solid_angle_sr: float) : None
        wrapper to update the location in all the populations
    estimate_field_population() : None
        wrapper to estimate the field output for each population
    do_kinematics(field_df: pandas.DataFrame) : pandas.DataFrame
        wrapper to call the kinematics generation in each population
    generate_fields() : None
        wrapper to generate all the populations
    write_astrotable(filename: str, df: pandas.DataFrame, extension: str) : None
        option to save the table as fits or VoTable.
    write_to_file(df: pandas.DataFrame) : str
        save the generated model, returns the filename
    get_filename(l_deg: float, b_deg: float, solid_angle_sr: float) : None
        generate the base of the filename (i.e. without extension) for a given location
    process_location(l_deg: float, b_deg: float,
            solid_angle_sr: float, save_data: bool) : Pandas.DataFrame Dict
        process a given location.
    process_all() : None
        process all locations as specified in the configuration
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the synthpop model.
        all inputs are forwarded to synthpop_utils.Parameters!
        """

        if not logger.file_logging_enabled:
            logger.setup_file_logging(20)

        # Initialize the model
        self.parms = sp_utils.Parameters(*args, **kwargs)
        # set random seed
        np.random.seed(self.parms.random_seed)

        # placeholder for population files
        self.population_files = None
        self.populations = []

        # placeholder for filename is updated every time a new location is set
        self.filename_base = None
        self.populations_are_initialized = False
        self.save_data = True

        if self.parms.post_processing_kwargs is None:
            self.post_processing = PostProcessing(self, logger).do_post_processing

        elif isinstance(self.parms.post_processing_kwargs, list):
            # have multiple post-processing class
            self.post_processing = []
            for post_kwargs in self.parms.post_processing_kwargs:
                post_kwargs.sun = self.parms.sun
                post_kwargs.model = self
                post_processing_class = sp_utils.get_subclass(
                    PostProcessing, post_kwargs)
                self.post_processing.append(post_processing_class.do_post_processing)

        else:
            # have single post-processing class
            self.parms.post_processing_kwargs.model = self
            post_processing_class = sp_utils.get_subclass(
                PostProcessing, self.parms.post_processing_kwargs)
            self.post_processing = post_processing_class.do_post_processing

        self.l_deg = None
        self.b_deg = None
        self.solid_angle = None
        self.solid_angle_unit = None

    def get_iter_loc(self) -> Iterator[Tuple[float, float]]:
        """ returns an iterator for the defined locations """
        return self.parms.location_generator()

    def init_populations(self, forced: bool = False) -> None:
        """
        Wrapper function that initializes all the populations for the model.
        Can initialize with coordinates, but it is not necessary.


        Parameters
        ----------
        forced : bool
            forced the initialisation, even if the populations are initialized before

        Returns
        -------

        """
        if self.populations_are_initialized & (not forced):
            # check if populations are already initialized
            # and preventing them from being initialized again
            logger.debug("populations are already initialized")
            return

        logger.create_info_section('initialize population')

        # get the model directory
        if not os.path.split(self.parms.model_name)[0]:
            models_dirname = os.path.join(self.parms.model_dir, self.parms.model_name)
            if not os.path.isdir(models_dirname):
                logger.critical(
                    f"{models_dirname} is not a directory")
                raise NotADirectoryError("Can not find the model directory for: "
                                         f"{self.parms.model_name}")
        elif os.path.isdir(self.parms.model_name):
            models_dirname = os.path.abspath(self.parms.model_name)
        else:
            logger.critical(f"{self.parms.model_name} is not a directory")
            raise NotADirectoryError("Can not find the model directory for: "
                                     f"{self.parms.model_name}")

        # make a list of all .json files in the model directory,
        # which will be all the included populations
        logger.info("read Population files from %s", self.parms.model_name)
        self.population_files = sorted(
            glob.glob(os.path.join(models_dirname, "*pop.json"))
            + glob.glob(os.path.join(models_dirname, "*popjson"))
            )

        if len(self.population_files) == 0:
            msg = f"No 'pop_json' file found in {models_dirname}, fall back to '.json' files"
            logger.warning(msg)
            self.population_files = sorted(
                glob.glob(os.path.join(models_dirname, "*.json"))
                )

        if len(self.population_files) == 0:
            msg = f"No Population file found in {models_dirname}"
            logger.critical(msg)
            raise FileNotFoundError(msg)

        # initialize the populations
        self.population_params = {pop_id: sp_utils.PopParams.parse_jsonfile(population_file)
            for pop_id, population_file in enumerate(self.population_files)
            }

        self.populations = [
            Population(pop_params, pop_id, self.parms)
            for pop_id, pop_params in self.population_params.items()
            ]
        logger.info("# all population are initialized")
        self.populations_are_initialized = True

    def update_location(
            self, l_deg: float, b_deg: float, solid_angle: float, solid_angle_unit: str = "deg^2",
            **positional_kwargs
            ) -> None:
        """
        Simple wrapper to update filename, logfile, and all
        populations in the model to the new location (l,b) and solid angle.

        Parameters
        ----------
        l_deg : float ['deg']
            galactic longitude for the center of the cone
        b_deg : float ['deg']
            galactic longitude for the center of the cone
        solid_angle : float [deg^2]
            steradians for the cone size
        solid_angle_unit : str
            unit for steradians
        positional_kwargs :
            any keyword to be passed to Position
        """

        if not self.save_data:
            self.filename_base = f'dump_file{np.random.randint(0, 999999):06d}'
        else:
            self.filename_base = self.get_filename(l_deg, b_deg, solid_angle)

        if not self.parms.overwrite:
            if os.path.isfile(ff := f"{self.filename_base}.{self.parms.output_file_type[0].lower()}"):
                msg = f"{ff} already exist!. Change the filename_pattern or use overwrite=True"
                logger.critical(msg)
                raise FileExistsError(msg)

        logger.update_location(f"{self.filename_base}.log", no_log_file=(not self.save_data))
        logger.create_info_section('update location')
        # update populations with the coordinates for the field
        logger.log(25, f"# set location to: ")
        logger.log(25, f"l, b = ({l_deg:.2f} deg, {b_deg:.2f} deg)")
        logger.log(25, f"# set solid_angle to:")
        logger.log(25, f"solid_angle = {solid_angle:.3e} {solid_angle_unit}")

        # placeholder for future position kwargs
        # (e.g. when we have the option for a pyramid cone)

        # update position in populations
        for population in self.populations:
            population.set_position(
                l_deg, b_deg, solid_angle, solid_angle_unit, **positional_kwargs)
        self.l_deg = l_deg
        self.b_deg = b_deg
        self.solid_angle = solid_angle
        self.solid_angle_unit = solid_angle_unit

        logger.debug("All Populations updated to (l,b) = "
                     f"({self.populations[0].b_deg:.2f},{self.populations[0].b_deg:.2f}) "
                     f"and Solid Angle = {self.populations[0].solid_angle_sr} sr.")

    def estimate_field_population(self) -> None:
        """
        A simple wrapper to estimate the total mass and number of stars to be drawn in a field
        """
        logger.debug('Estimate number of expected stars')
        for population in self.populations:
            mass_est, stars_est = population.estimate_field()
            logger.debug(f"Population {population.name} estimates: {mass_est:0.2f} M_sun "
                         f"will produce {stars_est:.1f} stars.")

    def do_kinematics(self, field_df: pandas.DataFrame) -> pandas.DataFrame:
        """

        A simple wrapper to call the kinematic generation of each population

        Note that the change is performed in place. I.e. It will overwrite
        what the kinematics in field_df is

        Parameters
        ----------
        field_df : DataFrame
           Generated stars

        Returns
        -------
        field_df
            dataframe of the generated stars
        """

        logger.info("# Determine velocities after all stars are evolved ")

        # extract columns
        all_mass = field_df.iloc[:, 4].to_numpy()
        all_dist, all_l_deg, all_b_deg = field_df.iloc[:, [6, 7, 8]].to_numpy().T
        all_x, all_y, all_z = field_df.iloc[:, [12, 13, 14]].to_numpy().T
        all_density_classes = tuple(pop.population_density for pop in self.populations)

        for pop_id, population in enumerate(self.populations):
            # select stars which belongs to population
            stars_with_pop_id = field_df[const.COL_NAMES[0]] == pop_id

            # call do_kinematics for the current population
            u, v, w, vr, mul, mub, vr_lsr = population.do_kinematics(
                dist=all_dist[stars_with_pop_id],
                star_l_deg=all_l_deg[stars_with_pop_id], star_b_deg=all_b_deg[stars_with_pop_id],
                x=all_x[stars_with_pop_id], y=all_y[stars_with_pop_id], z=all_z[stars_with_pop_id],
                mass=all_mass[stars_with_pop_id],
                all_x=all_x, all_y=all_y, all_z=all_z, all_mass=all_mass,
                all_density_classes=all_density_classes, pop_id=pop_id)

            # update velocities in population file
            field_df.loc[stars_with_pop_id, const.COL_NAMES[9]] = vr
            field_df.loc[stars_with_pop_id, const.COL_NAMES[10]] = mul
            field_df.loc[stars_with_pop_id, const.COL_NAMES[11]] = mub

            field_df.loc[stars_with_pop_id, const.COL_NAMES[15]] = u
            field_df.loc[stars_with_pop_id, const.COL_NAMES[16]] = v
            field_df.loc[stars_with_pop_id, const.COL_NAMES[17]] = w
            field_df.loc[stars_with_pop_id, const.COL_NAMES[18]] = vr_lsr

        return field_df

    def generate_fields(self) -> Tuple[pandas.DataFrame, Dict]:
        """
        calls the generate_field for all the populations and collects the data in a common
        pandas dataframe.

        Returns
        -------
        field_df: DataFrame
            DataFrame including the generated Stars from all populations
        distributions: dict
            collection of intrinsic distributions for each population
        """

        logger.create_info_section('Generate Field')

        # estimate number of stars, total-mass produced, etc
        self.estimate_field_population()
        # Placeholder to collect all the data_frames from the populations
        field_list = []
        distributions = {}

        for population in self.populations:
            # for each population, generate the field
            population_df, pop_distributions = population.generate_field()

            # collect distribution under the population name
            # be careful if two populations share the same name
            distributions[population.name] = pop_distributions

            logger.debug(
                "%s : Number of stars generated: %i (%i columns)",
                population.name, *population_df.shape)

            # collect data frame into field_list
            field_list.append(population_df)

        # combine them into one common data frame
        logger.create_info_section('Combine Populations')
        field_df = pandas.concat(field_list, ignore_index=True)

        logger.info('Number of stars generated: %i (%i columns)', *field_df.shape)
        # check if velocities should be generated after all positions are generated
        if self.parms.kinematics_at_the_end:
            field_df = self.do_kinematics(field_df)

        # check if faint stars and stars outside the grid should be kept or removed
        if self.parms.maglim[-1] != 'keep':
            logger.info('remove stars which are outside of the isochrone grid ')
            field_df = field_df[field_df[self.parms.maglim[0]] < self.parms.maglim[1]]
            logger.info('cleand field: Number of stars generated: %i (%i columns)', *field_df.shape)

        # log output columns and statistics
        logger.info(f"included_columns = {list(field_df.columns)}")
        sp_utils.log_basic_statistics(field_df, "stats_field")

        return field_df, distributions

    def write_astrotable(self, filename: str, df: pandas.DataFrame, extension: str) -> None:
        """
        save the result as fits file or votable by converting it to an astropy.table.Table object.

        Parameters
        ----------
        filename : str
            save location
        df : DataFrame
            Data frame of the generated model
        extension : str
            file extension. either "fits" or "vot"/"votable"

        Returns
        -------

        """

        # check if astropy was imported
        if 'astrotable' not in globals():
            logger.error("astropy is not installed, save as pickle instead")
            filename = f"{self.filename_base}.SAVE_ERROR.pkl"
            df.to_pickle(filename)
            return

        # convert vot to votable
        if extension == 'vot':
            extension = 'votable'

        # transfer to astropy table
        tab = astrotable.Table.from_pandas(df)
        # save table
        tab.write(filename, format=extension, overwrite=True)

    def write_to_file(self, df: pandas.DataFrame) -> str:
        """
        write the results to disc

        Parameters
        ----------
        df : DataFrame
            Model output

        Returns
        -------
        filename : str
            path of the created file save location.
        """
        # available formats:
        output_formatter = {
            'csv': [df.to_csv, 'csv', {'index': None, 'header': True, 'float_format': "%0.7e"}],
            'ssv': [df.to_csv, 'csv', {'index': None, 'header': True, 'float_format': "%0.7e", 'sep':' '}],
            'json': [df.to_json, 'json', {}],
            'html': [df.to_html, 'html', {'float_format': "%0.7e"}],
            'xml': [df.to_xml, 'xml', {}],
            'excel': [df.to_excel, 'xlsx', {'merge_cells': False, 'engine': None}],
            'hdf5': [df.to_hdf, 'h5', {"key": "data"}],
            'feather': [df.to_feather, 'ftr', {}],
            'parquet': [df.to_parquet, 'parquet', {}],
            'stata': [df.to_stata, 'dta', {}],
            'pickle': [df.to_pickle, 'pkl', {}],
            'sql': [df.to_sql, 'sql', {}],
            'fits': [self.write_astrotable, "fits", {"df": df, "extension": "fits"}],
            'votable': [self.write_astrotable, "votable", {"df": df, "extension": "votable"}],
            'vot': [self.write_astrotable, "vot", {"df": df, "extension": "vot"}],
            'ssv': [df.to_csv, 'csv', {'sep': ' ','index': None, 'header': True, 'float_format': "%0.7e"}],
            }
        output_file_type, output_save_kwargs = self.parms.output_file_type[:2]
        # get specific function for the defined format
        write_func, extension, kwargs = output_formatter.get(
            output_file_type.lower(), [df.to_pickle, 'SAVE_ERROR.pkl', {}]
            )
        kwargs.update(output_save_kwargs)
        # add extension to filename_base.
        filename = f"{self.filename_base}.{extension}"
        logger.log(25, 'write result to "%s"', filename)

        if extension.startswith('SAVE_ERROR'):
            logger.error(f"Invalid output file type {output_file_type}; "
                         f"save as pickle instead")

        # write dataframe to disk
        write_func(filename, **kwargs)

        return filename

    def get_filename(self, l_deg: float, b_deg: float, solid_angle: float) -> str:
        """
        create a file name for a given position

        Parameters
        ----------
        l_deg, b_deg : float [deg]
            galactic coordinates
        solid_angle: float [rad]
            solid angle of the cone

        Returns
        -------
        filename : str
            filename including location
        """
        if self.parms.scale_factor != 1:
            scale_factor_ending = f"_scaled{self.parms.scale_factor:.3f}"
        else:
            scale_factor_ending = ""

        file_keys = {
            "time": datetime.datetime.now().time(),
            "date": datetime.datetime.now().date(),
            "l_deg": l_deg,
            "b_deg": b_deg,
            "solid_angle_sr": solid_angle,
            "model_name": self.parms.model_name,
            "name_for_output": self.parms.name_for_output,
            }

        filename_base = os.path.join(
            self.parms.output_location,
            self.parms.output_filename_pattern.format(**file_keys) + scale_factor_ending)

        return filename_base

    def process_location(
            self, l_deg: float, b_deg: float, solid_angle: float, solid_angle_unit: str = 'deg^2',
            save_data: bool = True
            ) -> Tuple[pandas.DataFrame, Dict]:
        """
        Performs the field generation for a given position.

        Parameters
        ----------
        l_deg, b_deg : float [deg]
            galactic coordinates
        solid_angle : float [deg^2]
            Area of the cone
        solid_angle_unit : str
            Unit of the provided solid_angle
        save_data : bool
            If True the DataFrame is saved to disk
            If False the DataFrame are only returned,

        Returns
        -------
        field_df : DataFrame
            Generated stars as Pandas Dataframe
        distributions: dict
            collection of intrinsic distributions for each population
        """

        # store boolean if data should be stored to disc
        self.save_data = save_data

        # check if model has at least one population
        if len(self.populations) == 0:
            raise AttributeError("No initialized Population found! "
                                 "You must run 'init_populations()' prior to 'process_location()'!")

        # Step 1: Set the location and cone size
        self.update_location(l_deg=l_deg, b_deg=b_deg,
                             solid_angle=solid_angle, solid_angle_unit=solid_angle_unit)

        # Step 2: Generate all the fields for each population
        ti = time.time()
        field_df, distributions = self.generate_fields()
        t1 = time.time() - ti

        # Step 3: Save the results
        ti = time.time()
        if isinstance(self.post_processing, list):
            # have multiple post-processing
            for post_processing in self.post_processing:
                field_df = post_processing(field_df)
        else:
            # have single postprocessing
            field_df = self.post_processing(field_df)

        if self.save_data:
            logger.create_info_subsection('Save result')
            self.write_to_file(field_df)
        t2 = time.time() - ti

        # log end statement
        logger.log(15, f"{solid_angle}, {l_deg:.3f}, {b_deg:.3f}, done")
        logger.debug("---------------------------------------------------------------")
        logger.debug(f"Took total time {time.process_time()}")
        logger.debug(f'generate field: {t1:.1f}s | save field {t2:.1f}s')
        logger.info("---------------------------------------------------------------\n")

        return field_df, distributions

    def process_all(self, forced=False) -> None:
        """
        Calls the Initialization of all populations and calls the process for all given locations.
        """

        # Initializing all the different moduls.
        if (not self.populations_are_initialized) | forced:
            self.init_populations(forced=forced)

        # create output location, if it does not exist
        os.makedirs(self.parms.output_location, exist_ok=True)

        # Go through each b and l combination
        for l_b_deg in self.parms.loc:
            self.process_location(
                *l_b_deg, solid_angle=self.parms.solid_angle,
                solid_angle_unit=self.parms.solid_angle_unit)


def main(configfile: str = None, **kwargs):
    """ Run Synthpop in the standard mode """
    if configfile is None:
        # read input
        args = sp_utils.parser()

        # setup logger
        logger.setup_file_logging(20 + 5 * (args.quiet - args.verbose))

        # initialise model
        mod = SynthPop(
            specific_config=args.specific_config,
            default_config=args.default_config,
            model_dir=args.model_dir,
            **kwargs
            )
    else:
        mod = SynthPop(specific_config=configfile, **kwargs)

    # run model
    mod.process_all()

    # end logging
    logger.remove_file_handler()
