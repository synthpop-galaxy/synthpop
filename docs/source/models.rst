Models
======

A model is defined by a collection of populations. 
Each population is further described by the following 5 modules:
1. Population Density
2. Initial Mass Function
3. Age Distribution
4. Metallicity Distribution
5. Kinematics

A model is implemented as a dictionary containing separate json files, one for each population.
To avoid ambiguity, a ".popjson" extension should be used for these files. 
Each files must define the following keywords::

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
