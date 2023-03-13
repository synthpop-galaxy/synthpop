""" run Synthpop in the default way """
import os
if not os.path.isfile(os.path.join(os.path.dirname(__file__),'constants.py')):
    print("You need to migrate the to migrate the interactive part ")
    print("Pleas specify a Directory")
    from .migrate_interactive_part import migrate
    migrate()
del os
try:
    from .synthpop_main import main
except ImportError:
    from synthpop_main import main

main()
