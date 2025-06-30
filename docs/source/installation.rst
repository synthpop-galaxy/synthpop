Installation
============

To install SynthPop, you have two options:

1. Clone this repository and install all the requirements.
2. Use the command below to install using pip::

    pip install git+https://github.com/synthpop-galaxy/synthpop.git

If using pip, before you run SynthPop for the first time, you'll need to run::
    
    python -m synthpop.migrate_interactive_part

and follow the prompts to setup the interactive parts of the code in a user-accessible location.

You can instead specify the directory directly by running:: 

    python -m synthpop.migrate_interactive_part path_to_directory

Afterwards, SynthPop is ready to be used. 

This process can be undone by running:: 

    python -m synthpop.undo_migrate_interactive_part <path_to_directory>

with the path to directory being optional.

