"""
This file contains the Parameters Class.
Its task is to handel the configuration of SynthPop from the config files a Keyword input .
And provide the different task in a Namespace Object.
"""
__all__ = ["Parameters", "parser", "PopParams", "ModuleKwargs"]
__author__ = ["J. Klüter", "M.J. Huston"]
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2023-03-02"

import os
import json
import itertools
import logging
from typing import Iterator, Tuple, Dict, Optional, Union, List
import argparse
import numpy as np
import pydantic

if pydantic.__version__.startswith("2"):
    from pydantic import BaseModel, model_validator
else:
    from .pydantic_1_compatibility import BaseModel, model_validator
from .json_loader import json_loader
from .synthpop_logging import logger
from .sun_info import SunInfo

try:
    from constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR,
                           DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)
except (ImportError, ValueError):
    from ..constants import (SYNTHPOP_DIR, DEFAULT_MODEL_DIR,
                             DEFAULT_CONFIG_FILE, DEFAULT_CONFIG_DIR)


class ModuleKwargs(BaseModel, extra="allow"):
    name: Optional[str] = ""
    filename: Optional[str] = ""

    @classmethod
    @model_validator(mode="before")
    def validate_filename(cls, values):
        assert values.get("filename", "") != "" or values.get("name", "") != "", \
            "Either filename or name need to be specified"
        return values

    @property
    def init_kwargs(self):
        return {k: v for k, v in self.model_extra.items() if k != "json_file"}


class ExtLawKwargs(ModuleKwargs):
    bands: Optional[List[str]] = []
    min_lambda: Optional[float] = None
    max_lambda: Optional[float] = None
    lambdas: List[float] = None


class PopParams(BaseModel, extra="allow"):
    name: str

    imf_func_kwargs: ModuleKwargs
    age_func_kwargs: ModuleKwargs
    metallicity_func_kwargs: ModuleKwargs
    population_density_kwargs: ModuleKwargs
    kinematics_func_kwargs: ModuleKwargs

    evolution_kwargs: Optional[Union[ModuleKwargs, List[ModuleKwargs]]] = None

    av_mass_corr: Optional[float] = None
    n_star_corr: Optional[float] = None

    @classmethod
    def parse_jsonfile(cls, filename):
        try:
            data = json_loader(filename)
            obj = cls.model_validate(data)
        except Exception as e:
            logging.error(f"Error occurred when reading specifications from {filename}")
            raise e
        return obj


class Parameters:
    """
    NameSpace object to store the configuration parameters
    defined by keyword arguments, and in the configuration files
    """

    def __init__(
            self, specific_config: str or None = None,
            default_config: str = DEFAULT_CONFIG_FILE,
            model_dir: str = DEFAULT_MODEL_DIR,
            **kwargs
            ):
        """
        load configuration from files and kwarg arguments.

        Parameters
        ----------
        specific_config: str or None
            file of the specific configuration
        default_config: str
            file to the default configuration
        kwargs: dict
            set of keyword control arguments
        """
        self.l_set = None
        self.l_set_type = None
        self.b_set = None
        self.b_set_type = None

        logger.create_info_section('Settings')
        self._categories = {}
        # read settings form config files a kwargs arguments
        config_dir = self.read_default_config(default_config)
        self.read_specific_config(specific_config, config_dir)
        self.read_kwargs_config(kwargs)

        # generate random seed if not none
        # this allows to repeat the generation process later
        if not self.random_seed:
            self.random_seed = np.random.randint(0, 2 ** 31 - 1)

        # log settings to file
        self.log_settings()

        # check if Settings are ok
        if not self.validate_input():
            msg = "Settings Validation failed!." \
                  " Please ensure that all mandatory parameters are set."
            logger.critical(msg)
            raise ValueError(msg)

        # transfer l, b into a a location generator.
        self.loc = self.location_generator()
        self.model_dir = model_dir
        # check if model_base_name should be added to the output_location
        if self.output_location is None or self.output_location == "":
            self.output_location = os.path.join(SYNTHPOP_DIR, "outputfiles" + os.sep)

        if self.output_location.endswith(os.sep) or self.output_location.endswith("/"):
            self.output_location = os.path.join(self.output_location, self.name_for_output)

        if getattr(self, 'skip_lowmass_stars', False) \
                and getattr(self, 'kinematics_at_the_end', False):
            raise ValueError("'skip_lowmass_stars' and 'kinematics_at_the_end' "
                             "can not be set to true simultaneously")
        #
        self.sun = SunInfo(**self.sun, **self.lsr)
        # convert to ModuleKwargs BaseModels
        if isinstance(self.evolution_class, list):
            self.evolution_class = [
                ModuleKwargs.parse_obj(ev) for ev in self.evolution_class]
        else:
            self.evolution_class = ModuleKwargs.parse_obj(self.evolution_class)

        self.extinction_map_kwargs = ModuleKwargs.parse_obj(self.extinction_map_kwargs)

        if isinstance(self.extinction_law_kwargs, list):
            self.extinction_law_kwargs = [
                ExtLawKwargs.parse_obj(ext_law) for ext_law in self.extinction_law_kwargs
                ]
        else:
            self.extinction_law_kwargs = ExtLawKwargs.parse_obj(self.extinction_law_kwargs)

        if isinstance(self.post_processing_kwargs, list):
            self.post_processing_kwargs = [
                ModuleKwargs.parse_obj(pp) for pp in self.post_processing_kwargs]
        elif self.post_processing_kwargs is not None:
            self.post_processing_kwargs = ModuleKwargs.parse_obj(self.post_processing_kwargs)

        if isinstance(self.output_file_type, str):
            self.output_file_type = [self.output_file_type, {}]

    def validate_input(self):
        """ checks if all Mandatory files are provided"""
        out = True
        for key in self._categories["MANDATORY"]:
            if self.__getattribute__(key) is None:
                logger.critical(f"MANDATORY setting '{key}' is not defined in the configuration")
                out = False
        return out

    def location_generator(self) -> Iterator[Tuple[float, float]]:
        """
        converts l_set and b_set into a location generator object
        as defined by the l/b_set_type
        """
        if (self.l_set is None) or (self.b_set is None) or (self.solid_angle is None):
            logger.critical("Location or solid_angle_sr are not  defined in the settings! "
                            "Can not run main() or process_all()")

            # create

            def no_location():
                logger.critical("Location or solid_angle_sr are not defined in the settings! "
                                "Can not run main() or process_all()")
                for _ in []:
                    yield 0, 0

            return no_location()

        if (self.l_set_type == "pairs") or (self.b_set_type == "pairs"):
            return zip(self.l_set, self.b_set)

        if self.l_set_type == 'range':
            self.l_set = np.arange(*self.l_set)
            self.l_set_type = 'list'

        if self.b_set_type == 'range':
            self.b_set = np.arange(*self.b_set)
            self.b_set_type = 'list'

        return itertools.product(self.l_set, self.b_set)

    def rest_loc_iterator(self):
        self.loc = self.location_generator()

    def log_settings(self):
        """
        log the settings as json formatted file
        """
        logger.create_info_subsection('copy the following to a config file'
                                      ' to redo this model generation', 20)
        json_object = json.dumps(
            {key: item for key, item in self.__dict__.items() if not key.startswith('_')},
            indent=4)
        logger.info(json_object)

    def read_default_config(self, default_config_file: str):
        """
        reads settings from keyword arguments

        Parameters
        ----------
        default_config_file: str
            file and path of the default configuration file

        Return
        ------
        config_dir: str
            directory of default_config_file.

        """
        if os.path.isfile(default_config_file):
            config_dir = os.path.abspath(os.path.dirname(default_config_file))+'/'
            add_dir=''
        elif os.path.isfile(os.path.join(DEFAULT_CONFIG_DIR, default_config_file)):
            config_dir = DEFAULT_CONFIG_DIR+'/'
            add_dir=config_dir
        else:
            raise FileNotFoundError(f'default_config_file "{default_config_file}" was not found')

        logger.info("# reading default parameters from")
        logger.info(f"default_config_file =  {add_dir}{default_config_file} ")
        default = json_loader(add_dir+default_config_file)
        for cat, item in default.items():
            if cat.startswith("_"):
                continue
            self.__dict__.update(item)
            self._categories[cat] = item.keys()
        return config_dir

    def read_specific_config(self, config_file: str or None, config_dir: str = "."):
        """
        reads settings from keyword arguments
        Parameters
        ----------
        config_file: str
            file and path of the  configuration file
        config_dir: str
            directory where config files can be found
        """

        if config_file is None:
            return
        if not os.path.isfile(config_file):
            config_file_2 = os.path.join(config_dir, config_file)
            if not os.path.isfile(config_file_2):
                msg = f" can not find {config_file = } or {config_file_2}!"
                logger.critical(msg)
                raise FileNotFoundError(msg)

            config_file = config_file_2
        logger.info("# read configuration from ")
        logger.info(f"{config_file = !r} ")
        specified = json_loader(config_file)

        for cat, items in self._categories.items():
            spec_dict = specified.get(cat, specified)
            for item in items:
                if item in spec_dict:
                    self.__dict__.update({item: spec_dict[item]})

    def read_kwargs_config(self, kwargs: dict):
        """
        reads settings from keyword arguments

        Parameters
        ----------
        kwargs: dict
            dictionary of specifications
        """

        # for cat, items in self._categories.items():
        #     kwarg_dict = kwargs.get(cat, kwargs)
        #     for item in items:
        #         if item in kwarg_dict:
        #               self.__dict__.update({item:kwarg_dict[item]})
        self.__dict__.update(kwargs)


def parser() -> argparse.Namespace:
    """ synthpop control and help interface """

    parser_ = argparse.ArgumentParser(description="Running the SynthPop Model Generator")

    # Set control parameters
    parser_.add_argument(
        'specific_config', nargs='?',
        help="Config file to set parameter specialized for each model generation")
    parser_.add_argument(
        '--specific_config', '-s', dest='s', metavar='specific_config',
        help="option to pass a specific config")
    parser_.add_argument(
        'default_config', nargs='?',
        default=DEFAULT_CONFIG_FILE,
        help="Default config file to set common parameters (default is _default.synthpop_conf)")
    parser_.add_argument(
        '--default_config', '-d', dest='d', metavar='default_config',
        help="option to pass a default config")
    parser_.add_argument(
        '--model_directory', '-m', dest="model_dir",
        default=DEFAULT_MODEL_DIR,
        help="Path to the model directory")

    group = parser_.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="count", default=0)
    group.add_argument("-q", "--quiet", action="count", default=0)

    args = parser_.parse_args()

    # check if config is given as flag parameter
    if args.s is not None:
        if args.specific_config is not None:
            raise TypeError('got multiple values for argument "specific_config"')
        args.specific_config = args.s
    delattr(args, 's')
    if args.d is not None:
        if args.default_config is not parser.get_default('default_config'):
            raise TypeError('got multiple values for argument "default_config"')
        args.default_config = args.d
    delattr(args, 'd')

    return args
