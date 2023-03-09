# synthpop
Synthpop is a object oriented modular python framework to generate Synthetic Population Models. 
It will generate a star catalog following the specified model and configuration.




  


##  Model & Populations 
A Model is defined by a collection of Populations.
Each population is further described by the following 5 modules. 
1) Population Density 
2) Initial Mass Function
3) Age Distribution 
4) Metallicity Distribution
5) Kinematic

A Model is implemented as a dictonary containing for each Population a 
json files. These must define the following keywords. 

    "name" : "name_of_the_population"
    "imf_func_kwargs" : ...
    "age_func_kwargs" : ...
    "metallicity_func_kwargs" : ...
    "kinematics_func_kwargs" : ...
    "population_density_kwargs" : ...

each of the kwargs items includes a sub dictionary specifying the module (see below )

## Modules
In total 10 Moduls are used within Synthpop. These are: 
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

Each fulfilling a single task and can be specified independently.
Only the Isochrone System needs to be compatible with the Isochrone Interpolator. 
The predefined modules can be found in the modules directory tree. 
These are subclasses of a dedicated parent class. 
The parent classes can be found in the python files idicated with a "_"
For the Age class this can be found in 
    
    /modules/age/_age.py

We recommend to specify each subclass in a separate file within this file structure. 

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
Implementing a new Module is equivalent to creating a new subclass 
of the corresponding parent class. 
However, the footprint of the functions must not be changed.  


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
  



## run synthpop
### run SynthPop as individual script
  To run the SynthPop in the default mode:
  ```
  python synthpop config_file 
  ```
  this process all locations as defined in the config_file 
  Have a look at ```python synthpop -h ``` for additional arguments
  all results are saved at the location defined in the configuration. 
  (default: synthpop/output-files)
  
### import SynthPop to other script 
  importing synthpop to another script allows a bit more flexibility.
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
