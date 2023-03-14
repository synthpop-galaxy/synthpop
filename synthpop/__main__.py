""" run Synthpop in the default way """
import os
this_dir = os.path.dirname(__file__)
try:
    from .constants import SYNTHPOP_DIR

except ImportError:
    from .migrate_interactive_part import migrate
    migrate()

else:
    try:
        from .synthpop_main import main
    except ImportError:
        from synthpop_main import main
    main()
