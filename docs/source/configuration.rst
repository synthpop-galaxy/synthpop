SynthPop Configuration
=======================

SynthPop is configured via json files with a ``.synthpop_conf`` extension. Each run must have a default configuration, and may optionally have a specific configuration which overrides any default values.

The configuration sets the following: 
* Sightlines
* Coordinate system definition
* Population generation settings
* Extinction
* Isochrones + interpolator
* Photometry settings
* Output contents + formatting

Specific and Default Configurations
-----------------------------------
We provide a ``_default_config.synthpop_conf`` file, which is the default file used if a custom one is not entered.
A specific configuration will override any parameters defined in the default.

The default values are defined in config_files/_default_config.json.
Within the specific config file, it is sufficient to specify only the items 
which differ from the default config file.  
We recommend to include all information you think are 
usefull to understand the meaning of your generated catalog.
For example the config file can look as follows::
    
    {
      "model_base_name":"my_generated_model",
      "l_set": [0, 2],
      "b_set": [-1],
      "solid_angle":5e-7,
      "model_name":"besancon_Robin2003",
      "evolution_class":{
        "name":"MIST", 
        "interpolator":"CharonInterpolator"
        },
      "extinction_map_kwargs": {"name":"Marshall"}
      "extinction_law_kwargs": {"name":"ODonnellCardelli"}
    }

Setting Your Configuration
---------------------------
To set a configuration for a SynthPop run via command line, use::

  python -m synthpop <config_file>

or::

  python -m synthpop --specific_config <config_file> --default_config <default_config_file>

The equivalent process in a script that imports synthpop is::

  import synthpop
  model = synthpop.SynthPop('config_file',**kwargs)
  model.init_population()
  model.process_all()

where ``default_config`` is an optional keyword argument.

Configuration file formatting
-----------------------------

Our ``_default_config.json`` is sorted into different categories, which are not required but may be used for any configuration. Note that any arguments starting with an '#' are ignored, so we use these for comments. The optional category headings are: MANDITORY, SIGHTLINES, SEED, COORDINATE_SYSTEM, POPULATION_GENERATION, EXTINCTION_MAP, ISOCHRONE_INTERPOLATION, PHOTOMETRIC_OUTPUTS, and OUTPUT.

MANDITORY
^^^^^^^^^
**model_name**: string defining the model directory containing the population files

**name_for_output**: string defining the default directory and base for the output files

SIGHTLINES
^^^^^^^^^^
**l_set**

**l_set_type**

**b_set**

**l_set_type**

**solid_angle**

**solid_angle_unit**

SEED
^^^^

COORDINATE_SYSTEMS
^^^^^^^^^^^^^^^^^^

POPULATION_GENERATION
^^^^^^^^^^^^^^^^^^^^^

EXTINCTION_MAP
^^^^^^^^^^^^^^

ISOCHRONE_INTERPOLATION
^^^^^^^^^^^^^^^^^^^^^^^

PHOTOMETRIC_OUTPUTS
^^^^^^^^^^^^^^^^^^^

OUTPUT
^^^^^^
