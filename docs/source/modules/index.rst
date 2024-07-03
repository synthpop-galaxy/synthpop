Modules
=======

.. note::
  More module documentation coming soon

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
For instance, the parent class for the Age module can be found in::

.../modules/age/_age.py 

A uniform distribution subclasses is then specified in::

.../modules/age/uniform.py 

We recommend to specified each subclass in a separate file within this directory structure.
However, different locations are acceptable by specifying the filename, including the path . 


Define the used module
^^^^^^^^^^^^^^^^^^^^^^^
The usage of a module is either defined by the configuration or by a population file.
The used module is selected based on a dictionary:: 

    "used_module_kwargs":{
        "name" : "name_of_the_subclass"
        "filename" : "path_to_the_python_file"
        ... (further keys for the initialization)
        }

the filename is optional if the file is with named 
similar to the name of the subclass and in the default location. 
similar means either identical, all lower case, 
or as a snake_case filename for CamelCase class name.

Implementing a new Module
^^^^^^^^^^^^^^^^^^^^^^^^^^
Synthpop provides a variety of predefined modules to create synthetic populations, 
We encourage users to develop their own custom submodules to fit their needs.

To create a custom submodule, users simply can and define their own subclass of the appropriate parent class. 
Idealy in create a new Python file in the appropriate directory within the modules tree, 
Users can refer to the existing modules as a guide on how to structure their own custom module.
