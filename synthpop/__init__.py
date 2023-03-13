""" make Synthpop available for import """
import os
if not os.path.isfile(os.path.join(os.path.dirname(__file__),'constants.py')):
    print("You need to migrate the to migrate the interactive part ")
    print("Pleas specify a Directory")
    from .migrate_interactive_part import migrate
    migrate()
del os

from . import constants
from . import constants
from . import modules
from .synthpop_main import *
from .synthpop_main import __doc__