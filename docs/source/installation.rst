Installation
============

To install SynthPop, you have two options:

1. Clone this repository and install all the requirements.
2. Use the command below to install using pip::

    pip install git+https://github.com/synthpop-galaxy/synthpop.git

When using SynthPop for the first time, you are ask to specify a directory to 
store all files you might want to interact with. 
These are The module, models, constants. It will also be the default location to store the isochrone data 
and output files. 

A simple GUI will be used specify the directory.

You can instead specify the directory directly by running:: 

    python -m synthpop.migrate_interactive_part path_to_directory

Afterwards, SynthPop is ready to be used. 
