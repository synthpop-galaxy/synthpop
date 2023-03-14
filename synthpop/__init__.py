""" make Synthpop available for import """
import sys

try:
    from .constants import SYNTHPOP_DIR

except ImportError:
    if sys.argv[0] != '-m':
        raise ImportError("You must migrate the files first using migrate_interactive_part")

else:
    from . import modules
    from . import constants
    from .synthpop_main import *
    from .synthpop_main import __doc__
