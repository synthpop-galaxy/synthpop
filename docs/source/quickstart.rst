Quickstart
=======================

Installation
------------
To install SynthPop, download the repository.
If you plan to import SynthPop to other scrips 
ensure that the parent dictionary is in your Python Path.

Requirements
^^^^^^^^^^^^
The following packages are needed

 * numpy
 * scipy
 * pandas
 * request
 * json
 * tqdm
 * tables
 * dustmaps
 * astropy
  
The tables, dustmaps, and astropy packages are optional but highly recommended.
Install all by using::

  pip install -r requirements.txt

Run SynthPop
------------
Synthpop can either be run as an individual script 
or can be imported to other scripts.

Run SynthPop as individual script
---------------------------------
To run the SynthPop type::
  
  python synthpop config_file 
  
This processes all locations as defined in the config_file.
Have a look at ``python synthpop -h `` for additional arguments.
All results are stored in the defined output location
(default: synthpop/output-files).
  
Import SynthPop to other script 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
SynthPop can be imported to another script, which provides more flexibility. 
to do so, ensure that the parent directory is within your python Path::
  
  import synthpop
  
  model = synthpop.SynthPop('config_file',**kwargs)
  model.init_population()
  
All attributes in the configuration can also be modified by a keyword argument.
  
To run all specified locations::
  
  model.process_all() 

To run a specified locations::
  
  data, distributions = model.process_location(
        l_deg, b_deg,solid_angle) 
  
while ``process_all()`` only saves the results to disk,
``process_location()`` also returns the dataframe and several distributions
(currently only distance distributions for each population)

Configure SynthPop
------------------
See :ref:`Configuration`
 

Select magnitudes 
^^^^^^^^^^^^^^^^^
The ``"chosen_band"`` argument is used to specify the returned filters, 
With the default MIST isochrone system, any available filter can be selected. 
However, these have to follow the MIST column names (e.g. The PanSTARRS g magnitude can only specify by "PS_g")
Note that the implement MIST module downloads the data automatically. 

