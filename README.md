# Synthpop

(For full documentation, see our ReadTheDocs)

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

A model is implemented as a dictionary containing separate json files, one for each population.
To avoid ambiguity, a ".popjson" extension should be used for these files. 
Each files must define the following keywords:


    "name" : "name_of_the_population"
    "imf_func_kwargs" : ...
    "age_func_kwargs" : ...
    "metallicity_func_kwargs" : ...
    "kinematics_func_kwargs" : ...
    "population_density_kwargs" : ...

Each of the kwargs items includes a sub-dictionary 
specifying the module (see below).
Additional kwargs can be used to define a population specific
evolutionary modul and galactic warp.

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

Each module fulfills a single task and can be independently specified.
Only the Isochrone System needs to be compatible with the Isochrone Interpolator.
All predefined modules are located in the modules directory tree 
and are subclasses of a dedicated parent class, which is denoted by a "_" in the name of the corresponding Python file.
For instance, the parent class for the Age module can be found in:\

``.../modules/age/_age.py ``

A uniform distribution subclasses is then specified in:

``.../modules/age/uniform.py ``

We recommend to specified each subclass in a separate file within this directory structure.
However, different locations are acceptable by specifying the filename, including the path . 


### Define the used module. 
The usage of a module is either defined by the configuration or by a population file.
The used module is selected based on a dictionary. 

    "used_module_kwargs":{
        "name" : "name_of_the_subclass"
        "filename" : "path_to_the_python_file"
        ... (further keys for the initialization)
        }

the filename is optional if the file is with named 
similar to the name of the subclass and in the default location. 
similar means either identical, all lower case, 
or as a snake_case filename for CamelCase class name.

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

To install synthpop, you have two options:

Clone this repository and install all the requirements.
Use the command below to install using pip:

```
pip install git+https://github.com/synthpop-galaxy/synthpop.git
``` 
When using SynthPop for the first time from a pip install, you will need to run the following command to move a few interactive files and directories to a user-accessible location:

```
python -m synthpop.migrate_interactive_part
```
Only afterwards, synthpop is ready to be used. 



## Use SynthPop
### Run Synthpop as individual script
To run SynthPop in the default mode, use the following command:
  ```
  python -m synthpop config_filename 
  ```
This processes all locations as defined in the config_filename. 
The config_filename should either be in the ``config_file`` directoryor should include the complete path.
As an example, you can use the predifined ``my_config.synthpop_conf`` file. 
Note that you do not need to include a ``-m`` flag when SynthPop is within your current working directory. For additional arguments, see python ``-m synthpop -h``.

The generated catalogs are saved at the location defined in the configuration (default: your_specified_directory/output-files).


  
### Import SynthPop to other script 
Importing SynthPop to another script allows for more flexibility. 
To do so, ensure that the parent directory is within your Python path and use the following code:
  ```
  import synthpop
  
  model = synthpop.SynthPop(config_file_name, **kwargs)
  model.init_population()
  ```
All attributes of the configuration can also be specified by a keyword argument. 
It is then possible to run all specified locations via
  ```
  model.process_all() 
  ```
or to run a specified location only:

  ```
  data, distributions = model.process_location(
        l_deg, b_deg, solid_angle, save_data=True) 
  ```
While ``process_all()`` only saves the results to disk, ``process_location()`` also returns the dataframe and several distributions (currently only distance distributions for each population). Later als have a ``save_data``- flag. If it is set to False, the data will only be returned and not saved.


## Acknowledge Synthpop 
  If you think SynthPop was usfull for you work, please cite Klüter et al. (in prep). 
  Please also include citations for key components of the generation process. 
  These includes, but is not limited, to the used model, isochrone system and extinction map.

## Getting in touch:
  If users encounter any issues while using Synthpop during the development/implementation of new submodules or models, 
  they can reach out to the development team through the GitHub issue tracker. 
  We welcome any feedback, bug reports, or feature requests that can help improve Synthpop.
