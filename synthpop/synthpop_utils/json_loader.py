"""
loader of json file,
which ignores  comments starting with "#" and allows nested files
whenever the "json_file" keyword is detected
"""

__all__ = ["json_loader", ]
__author__ = "S. Johnson"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-07-06"

import os
import sys
import json
import logging

try:
    import constants as const
except (ImportError, ValueError):
    from .. import constants as const

logger = logging.getLogger('synthpop')


def scrub_dict(obj, bad_key="_this_is_bad"):
    """
    Method to remove comments from a dictionary, where comments are keys that 
    have a hash symbol ('#') anywhere in them.

    Adapted from https://stackoverflow.com/a/20692955/17159462
    """
    if isinstance(obj, dict):
        # the call to `list` is useless for py2 but makes
        # the code py2/py3 compatible
        for key in list(obj.keys()):
            if bad_key in key:
                del obj[key]
            else:
                scrub_dict(obj[key], bad_key)
    elif isinstance(obj, list):
        for i in reversed(range(len(obj))):
            if obj[i] == bad_key:
                del obj[i]
            else:
                scrub_dict(obj[i], bad_key)

    else:
        # neither a dict nor a list, do nothing
        pass


def substitutes_files(data, filename):
    """
    searches for json_file keywords in the dictionary, and loads the connected files.

    Parameters
    ----------
    data: dict
        dictionary where files should be substituted
    filename
        name of the original loading sequence
        used to track the location

    Returns
    -------

    """
    json_file = data.get('json_file')
    # if json_file is None:
    for key, value in data.items():
        if isinstance(value, dict):
            substitutes_files(value, filename)
    if json_file is None:
        return
        # check different lo
    for directory in ['', os.path.dirname(filename), const.SYNTHPOP_DIR]:
        file = os.path.join(directory, json_file)
        if os.path.isfile(file):
            break
    else:
        msg = f"{json_file} can not found (specified in {filename})"
        logger.critical(msg)
        raise FileNotFoundError(msg)

    with open(file) as handle:
        data.update(json.load(handle))


def json_loader(filename, json_file_branch=None):
    """
    Custom JSON loader that removes comments, substitutes filenames.

    Parameters
    ----------
    filename : str
        filename of the json file implemented
    json_file_branch : list
        used to detect recursions
    Returns
    -------
    data: dictionary
        loaded json file with removes comments, substitutes filenames.

    """
    if json_file_branch is None:
        json_file_branch = []

    if filename in json_file_branch:
        msg = f"infinite loop detected in:\n {' -> '.join(json_file_branch)} -> {filename}"
        logger.critical(msg)
        raise RecursionError(msg)

    # open up the file
    with open(filename) as handle:
        # format the json into data
        data = json.load(handle)

    data['_filename'] = filename

    # check to see if any of the top level keys indicate a JSON file for kwargs
    # only key should be "json_file"
    substitutes_files(data, filename)

    # recursively remove any keys with '#'
    # this operates on the dict itself, so nothing to return
    scrub_dict(data, '#')

    return data
