# Synthpop

Synthpop is an object-oriented, modular Python framework 
for generating synthetic population models. 
It generates a star catalog following the specified model and configuration.

## Model & Populations

A model is defined by a collection of populations. 
Each population is further described by the following 5 modules:
1. Population Density
2. Initial Mass Function
3. Age Distribution
4. Metallicity Distribution
5. Kinematics

A model is implemented as a dictionary containing, 
for each population, a JSON file. 
These files must define the following keywords:


    "name" : "name_of_the_population"
    "imf_func_kwargs" : ...
    "age_func_kwargs" : ...
    "metallicity_func_kwargs" : ...
    "kinematics_func_kwargs" : ...
    "population_density_kwargs" : ...

Each of the kwargs items includes a sub-dictionary 
specifying the module (see below).

## Modules

Synthpop employs 10 modules to perform various tasks, namely:

1) Population Density
2) Initial Mass Function
3) Age Distribution
4) Metallicity Distribution
5) Kinematic
6) Isochrone System
7) Isochrone Interpolator
8) Extinction Map
9) Extinction Law
10) Post-Processing

11) Each module fulfills a single task and can be independently specified.
Only the Isochrone System needs to be compatible with the Isochrone Interpolator.
All predefined modules are located in the modules directory tree 
and are subclasses of a dedicated parent class, which is denoted by a "_" in the name of the corresponding Python file.
For instance, the parent class for the Age module can be found in:

``.../modules/age/_age.py ``

We recommend that each subclass be specified in a separate file within this directory structure.

### Define the used Module. 
The usage of a module is either defined by the configuration or by a population file.
The used module is selected based on a dictionary. 

    "used_module":{
        "name" : "name_of_the_subclass"
        "filename" : "path_to_the_python_file"
        ... (further keys for the initialization)
        }

the file_name is optional if the file is with named 
similar to the name of the subclass and in the default location. 
similar means either identical, all lower case, 
or as snake_case for  CamelCase class name.  

### Implementing a new Module
Synthpop provides a variety of predefined modules to create synthetic populations, 
We encourage users to develop their own custom submodules to fit their needs.

To create a custom submodule, users simply can and define their own subclass of the appropriate parent class. 
Idealy in create a new Python file in the appropriate directory within the modules tree, 
Users can refer to the existing modules as a guide on how to structure their own custom module.



## Configuration 
  SynthPop is controlled by a config json file. 
  Within these you can specify the used model, 
  the isochrone system & interpolator extinction
  the wanted Properties for the output, output location 
  and several other control keys.
  The default values can be found in config_files/_default_config.json.
  Within the config file it is sufficient to specify only the items 
  which differs from the default_config file.  
  We recommend to include all information you think are 
  usefull to understand the meaning of your generated catalog.
  For example the config file can look as follows:
    
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
  

## Installation
To install synthpop, you can either clone this repository, and install all the requirements. 
or you can use 
```
pip install git+https://github.com/synthpop-galaxy/synthpop.git 
``` 
In this case you might want to migrade the model, module and conifgurations to a easyly accessable directory. 
To do so, you need to run 
```
python -m synthpop.migrate_interactive_part {path_to_directory}
```
You can also used the buildin gui to select a directory. 


## Use SynthPop
### Run Synthpop as individual script
  To run the SynthPop in the default mode:
  ```
  python synthpop config_filename 
  ```
  this process all locations as defined in the config_filename. 
  The config_filename sould either be in the ``config_file`` directory. 
  Or should include the complete path.
  Note you need to include a "-m" flag when installed by pip install.
  For additional arguments see `python synthpop -h `. 

  The generated catalogs are saved at the location defined in the configuration. 
  (default: synthpop/output-files)
  
### Omport SynthPop to other script 
  Importing synthpop to another script allows a bit more flexibility.
  To do so, ensure that the parent directory is within your python path, and use 
  ```
  import synthpop
  
  model = synthpop.SynthPop(config_file_name, **kwargs)
  model.init_population()
  ```
  all attributes of the configuration can also specified by a keyword argument
  It is then possible  to run all specified locations via 
  ```
  model.process_all() 
  ```
  or to run a specified location only: 
  ```
  data, distributions = model.process_location(
        l_deg, b_deg, solid_angle, save_data=True) 
  ```
  while ```process_all()``` only saves the results to disk,
  ```process_location()``` also returns the dataframe and several distributions
  (currently only distance distributions for each population)
  if save_data is False, the data will only be returned and not be saved. 


## Acknowledge Synthpop 
  If you think SynthPop was usfull for you work, 
  Please cite Klüter et al. (in prep), and include citations for key components of the generation process. 
  These includes, but is not limited, to the used model, isochrone system and extinction map.

## Getting in touch:
  If users encounter any issues during the development of their custom submodules or model or while using Synthpop, 
  they can reach out to the development team through the GitHub issue tracker. 
  We welcome any feedback, bug reports, or feature requests that can help improve Synthpop.