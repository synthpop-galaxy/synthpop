""" make Synthpop available for import """
try:
    from . import constants
except ImportError:
    from . import migrate_interactive_part
else:
    from . import modules
    from .synthpop_main import *
    from .synthpop_main import __doc__