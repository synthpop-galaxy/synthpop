Installation
============

To install synthpop, you have two options:

1. Clone this repository and install all the requirements.
2. Use the command below to install using pip::

    pip install git+https://github.com/synthpop-galaxy/synthpop.git

When using SynthPop for the first time, you are ask to specify a directory to 
store all files you might want to interact with. 
These are The module, models, constants. It will also be the default location to store the isochrone data 
and output files. 

A simple GUI will be used specify the directory
You can also specify the directory directly
To do so, run the following command:: 

    python -m synthpop.migrate_interactive_part path_to_directory

Only afterwards, synthpop is ready to be used. 
