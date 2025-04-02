"""
This file included ths SubClass selection and importing process 
"""
__all__ = ["get_subclass"]
__author__ = "Jonas Klüter"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2023-03-01"

import sys
import os
import importlib
import re
import inspect
from typing import Type, Dict, TypeVar, Union, Callable, List
from .synthpop_logging import logger
from .synthpop_control import ModuleKwargs
T = TypeVar('T')


class SubClassLoader:
    """
    This class returns a callable that handles
    the subclass detection loading and initialization process
    It takes track of previously imported moduls.
    """

    def __init__(self):
        # Location to store modules . This prevents loading the same file multiple times.
        self._loaded_submodules = {}
        self.do_logging = True

    def __call__(
            self,
            ParentClass: Type[T],
            modul_kwargs: ModuleKwargs,
            initialize: bool = True,
            keyword_name: str = 'kwargs',
            population_file: str = '??',
            no_logging: bool = False
            ) -> Union[Type[T], T]:
        """
            A method to get a subclass of ParentClass with either the name or filename
            that contains the subclass. Then, pass **kwargs to create an instance of
            the subclass.

            Parameters
            ----------
            ParentClass : class
            kwargs : dict
                keywords for the initialisation must have either a
                name or filename keyword to used to specify the subclass
                not_a_sub_class can be used to specify that a class is not written as subclass
                it must have the same functionality.

            initialize : bool, optional
                return initialized instance
            keyword_name : str
                string for a better logging information
            population_file: str
                origin for of the kwargs for better logging information
            no_logging : bool
                do not log information
            Returns
            -------
            subcls : class
                Found subclass, depending on the initialize keyword it is initialized or not.
            """

        # update logging information
        self.do_logging = not no_logging
        # get subclass_name, filename and not_a_sub_class  from kwargs
        # if isinstance(modul_kwargs, dict):
        #     modul_kwargs = ModuleKwargs.parse_obj(modul_kwargs)
        kwargs = modul_kwargs.init_kwargs
        subclass_name = modul_kwargs.name
        filename = modul_kwargs.filename
        not_a_sub_class = kwargs.pop('not_a_sub_class', False)

        if (not subclass_name) and (not filename):
            msg = f"neither a name or filename is given (from {population_file})"
            logger.critical(msg)
            raise AttributeError(msg)

        # estimate default location for ParentClass
        parent_dir = os.path.relpath(os.path.dirname(inspect.getfile(ParentClass)))

        if filename == '':
            # if determine filename from name, if not defined.
            filename = self.get_filename(subclass_name, parent_dir, population_file)
        else:
            # check if file exist
            # also add path if it is  not in included in  filename
            filename = self.check_file_location(filename, parent_dir, population_file)

        # load filename as submodule
        submodule = self.get_submodule(filename, ParentClass.__module__)

        # detect subclasses in that file
        subclasses = self.get_subclasses_from_module(
            submodule, subclass_name, ParentClass, not_a_sub_class
            )

        # check if a unique subclass was detected
        subcls = self.check_unambiguity(
            subclasses, ParentClass, subclass_name, filename, population_file
            )

        if initialize:
            # initialize subclass
            return self.initialize_class(subcls, kwargs, keyword_name)
        else:
            return subcls

    def get_filename(self, name: str, parent_dir: str, population_file) -> str:
        """
        estimates the filename from the class name
        it checks three different conventions. For a given name: 'ThisIsAClass' these are:
        1 same as class name: ThisIsAClass.py
        2 all lower letters: thisisaclass.py
        3 convert CamelCase to snake_case: this_is_a_class.py
        They are prioritized in the given order.
        """

        if check_file_case_sensitive(filename := os.path.join(parent_dir, name + '.py')):
            if self.do_logging:
                logger.debug('found file: %s', filename)
            return filename

        if check_file_case_sensitive(filename := os.path.join(parent_dir, name.lower() + '.py')):
            if self.do_logging:
                logger.debug('found file: %s', filename)
            return filename

        # translate CamelCase to snake_case
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        snakecase = re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
        if check_file_case_sensitive(filename := os.path.join(parent_dir, snakecase + '.py')):
            if self.do_logging:
                logger.debug('found file: %s', filename)
            return filename

        msg = f'No subclasses found with name = {name} (from {population_file})'
        logger.critical(msg)
        raise NameError(msg)

    def check_file_location(self, filename, parent_dir, population_file):
        """
        check if a file exist.
        find the file if no path is specified
        it can be either at the default location or in the current directory.
        """
        if os.path.dirname(filename) != '':
            if not check_file_case_sensitive(filename):
                msg = "can not find {filename}"
                logger.critical(msg)
                raise FileNotFoundError(filename)
            return filename

        if check_file_case_sensitive(os.path.join(parent_dir, filename)):
            if self.do_logging:
                logger.debug(f"{filename} was found at the default location {parent_dir}")
            return os.path.join(parent_dir, filename)

        if check_file_case_sensitive(filename):
            if self.do_logging:
                logger.debug(f"{filename} was found in the current directory")
            return filename

        msg = f"No such file {filename} (from {population_file})"
        logger.critical(msg)
        raise FileNotFoundError(msg)

    def get_submodule(self, filename, module):
        """
        load the module from a given file
        stores them in a dictionary for multiple uses
        """

        if filename in self._loaded_submodules:
            submodule = self._loaded_submodules.get(filename)
        else:
            spec = importlib.util.spec_from_file_location(module, filename)
            submodule = importlib.util.module_from_spec(spec)
            try:
                spec.loader.exec_module(submodule)
            except ModuleNotFoundError:
                sys.path.insert(0, os.path.dirname(filename))
                spec.loader.exec_module(submodule)
                sys.path.pop(0)

            self._loaded_submodules.update({filename: submodule})

        return submodule

    def get_subclasses_from_module(self, submodule, name, ParentClass, not_a_sub_class=False):
        """

        Parameters
        ----------
        submodule : loaded module from file
        name : str
            name of the subclass, can be ''
        ParentClass : class

        not_a_sub_class : bool
            if True it will not be checked if the class is a subclass

        Returns
        -------
        subcls : list[class]
            list of detected subclasses
        """

        # allow class to not be a subclass
        if not_a_sub_class:
            if self.do_logging:
                logger.debug('Do not check if class is a subclass')
            # check if attribute exist and is a class
            if hasattr(submodule, name):
                obj = getattr(submodule, name)
                if isinstance(obj, type):
                    # return as list to have the same footprint.
                    return [obj]
            # raise error if class was not found
            msg = f"submodule {submodule.__file__} has no object attribute named {name}"
            logger.critical(msg)
            raise AttributeError(msg)

        name_camel = ''.join(ele.title() for ele in name.split('_'))

        subcls = []
        subcls_weak = []

        for obj_name in dir(submodule):
            if (name != '') and (name != obj_name) and (name_camel != obj_name):
                continue
            obj = getattr(submodule, obj_name)

            if isinstance(obj, type):
                if obj in ParentClass.__subclasses__():
                    if name == obj_name:
                        subcls.append(obj)
                    else:
                        subcls_weak.append(obj)
                elif obj.__base__.__name__ == ParentClass.__name__:
                    # check if all functions in ParentClass are also defined in obj
                    if all(
                            hasattr(obj.__base__, func) for func in dir(ParentClass) if
                            not func.startswith("__")
                            ):
                        subcls_weak.append(obj)

        # do some logging
        if len(subcls) != 0:
            return subcls
        else:
            return subcls_weak

    def check_unambiguity(
            self,
            subclasses: List[Type[T]], ParentClass: Type[T],
            name: str = '??', filename: str = '??', population_file: str = '??'
            ) -> Type[T]:
        """
        check if the exactly 1 subclass was found
        Parameters
        ----------
        subclasses : list[class]
            list of subclasses
        ParentClass : class
        name : str
            name of the subclass for better logging messages
        filename : str
            filename of the module for better logging messages
        population_file : str
            population file that included the kwargs. Used for a better logging messages
        Returns
        -------
        subclass: class
            single subclass
        """
        if len(subclasses) > 1:
            msg = f"Found more than one subclasses with name = '{name}',  filename = '{filename}'" \
                  f" (from {population_file})"
            logger.critical(msg)
            logger.critical(subclasses)
            raise ModuleNotFoundError(msg)

        if len(subclasses) == 0:
            msg = f"No subclasses found with name = '{name}', filename = '{filename}'" \
                  f" (from {population_file})"
            logger.critical(msg)
            logger.critical(f"'{ParentClass.__subclasses__()}'")
            raise ModuleNotFoundError(msg)

        if self.do_logging:
            logger.log(15, f"set {ParentClass.__name__} class to {subclasses[0].__name__}")

        return subclasses[0]

    def initialize_class(self, sub_class, kwargs, keyword_name):
        """ initialize class with keywords """
        if self.do_logging:
            logger.debug(f'initialize class with keywords {kwargs}')
            msg = f'"{keyword_name}" : {kwargs}'.replace("'", '"')
            msg = msg.replace('True', 'true').replace('False', 'false')
            msg = msg.replace('None', 'null')
            logger.log(15, msg)
        return sub_class(logger=logger, **kwargs)


def check_file_case_sensitive(file: str) -> bool:
    """
    This function forces a case-sensitive check even on insensitive filesystems
    Parameters
    ----------
    file: str

    Returns
    -------

    """
    if not os.path.isfile(file):
        return False
    directory, filename = os.path.split(file)
    return filename in os.listdir(directory)


get_subclass = SubClassLoader()
