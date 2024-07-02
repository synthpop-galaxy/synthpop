Configuration
=======================

.. note::
    Page under construction

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
**model_name**: (string) model directory containing the population files

**name_for_output**: (string) default directory and base for the output files

SIGHTLINES
^^^^^^^^^^
**l_set**, **b_set**

**l_set_type**, **b_set_type**

**solid_angle**

**solid_angle_unit**

SEED
^^^^
**random_seed**

COORDINATE_SYSTEMS
^^^^^^^^^^^^^^^^^^
**sun**: dictionary containing the following values:
* **x**, **y**, **z**: location of the Sun in cartesian Galactic coordinates (kpc)
* **u**, **v**, **w**: motion of the Sun in cartesian Galactic coordinates (km/s)
* **l_apex_deg**, **b_apex_deg**: direction of the Solar apex in Galactic coordinates (degree)

**lsr**: dictionary containing the following values:
* **u_lsr**, **v_lsr**, **w_lsr**: velocity of the local standard of rest (km/s) in cartesian Galactic coordinates

**warp**: dictionary describing the warp of the galaxy (can be overwritten in population files)
* **r_warp**: radius where the warp begins (kpc, R_W in equation below)
* **amp_warp**: warp amplitude (a in equation below) [note: set this to null to use **amp_warp_pos**+**amp_warp_neg**]
* **amp_warp_pos**, **amp_warp_neg**: warp amplitude (a in equation below), for positive and negative phi [note: set these to null to use **amp_warp**]
* **alpha_warp**: exponent of warp power law (b in equation below)
* **phi_warp_deg**, **phi_warp_rad**: angle for line of notes, where one should be specified with indicated unit (degree or radian)
where warp is parametrized as z_w=a(R-R_W)^b sin(phi-phi_W) (see `Chen et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019NatAs...3..320C/abstract>`_)

POPULATION_GENERATION
^^^^^^^^^^^^^^^^^^^^^
**max_distance**

**distance_step_size**

**window_type**: dictionary containing the following:
* **window_type**
* **kwargs**

**mass_lims**: range of initial stellar masses to produce

**N_mc_totmass**: number of stars to use to estimate number of stars needed per slice

**

EXTINCTION_MAP
^^^^^^^^^^^^^^

ISOCHRONE_INTERPOLATION
^^^^^^^^^^^^^^^^^^^^^^^

PHOTOMETRIC_OUTPUTS
^^^^^^^^^^^^^^^^^^^

OUTPUT
^^^^^^
