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

We provide a ``_default_config.synthpop_conf`` file, which is the automatic default file used if a custom one is not entered. 

Setting Your Configuration
---------------------------
To set a configuration for a SynthPop run via command line, use::

  python synthpop config_file



The default values are defined in config_files/_default_config.json.
Within the config file, it is sufficient to specify only the items 
which differ from the _default_config file.  
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

Note that _default_config.json is sorted into different categories, 
These can be (but don't need to be) translated into the config file. 
Also note that arguments starting with an '#' are ignored.


