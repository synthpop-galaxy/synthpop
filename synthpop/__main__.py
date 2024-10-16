
""" run Synthpop in the default way """
flag = True
import_with_dot = True
try:
    from .constants import SYNTHPOP_DIR

except ImportError as e:
    if str(e) == 'attempted relative import with no known parent package':
        import_with_dot = False
        try:
            from constants import SYNTHPOP_DIR
        except ImportError:
            flag = False
    else:
        flag = False

if flag:
    if import_with_dot:
        from .synthpop_main import main
    else:
        from synthpop_main import main
    main()

else:
    if import_with_dot:
        from .migrate_interactive_part import migrate
    else:
        from migrate_interactive_part import migrate

    migrate()
