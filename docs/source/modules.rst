Modules
=======

.. note::
  More module documentation coming soon

Synthpop employs 10 modules to perform various tasks:

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

A uniform distribution subclass is then specified in::

.../modules/age/uniform.py 

We recommend to specify each subclass in a separate file within this directory structure.
However, different locations are acceptable by specifying the filename, including the path when referring to the modules in your configuration and population json files. 

Selecting Modules
-----------------
The usage of a module is either defined by the configuration or by a population file.
The used module is selected based on a dictionary:: 

    "<module_type>_kwargs":{
        "name" : "name_of_the_subclass",
        "filename" : "path_to_the_python_file",
        "<additional_kwargs>" : <value>
        }

The filename is optional if the file is in the default location and named similarly, i.e. identical, same but all lower case, 
or as a snake_case filename for CamelCase class name.

Implementing New Modules
------------------------
Synthpop provides a variety of predefined modules to create synthetic populations. 
We encourage users to develop their own custom submodules to fit their needs.

To create a custom submodule, users simply can and define their own subclass of the appropriate parent class. 
This should be done in a new Python file in the appropriate directory within the modules tree, referring to existing modules for structure.

Well-tested modules may be added to the SynthPop repository. Alternatively, users are welcome to share modules and other SynthPop-related code and data via `Zenodo <https://zenodo.org/communities/synthpop/records?q=&l=list&p=1&s=10&sort=newest>`_, which we encourage for any publications using SynthPop.

Additional Module Information
-----------------------------
See the synthpop.modules API autodocumentation for more details on individual modules.
