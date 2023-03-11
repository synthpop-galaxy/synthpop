""" run Synthpop in the default way """
try:
    from .synthpop_main import main
except ImportError:
    from .synthpop_main import main


main()
