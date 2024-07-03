Usage
============
  
Run Synthpop as individual script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To run SynthPop in the default mode, use the following command::

    python -m synthpop config_filename 

This processes all locations as defined in the config_filename. 
The config_filename should either be in the ``synthpop/config_file`` directory or should include the complete path.
As an example, you can use the predifined ``my_config.synthpop_conf`` file. 
Note that you do not need to include a ``-m`` flag when SynthPop is within your current working directory. For additional arguments, see ``python -m synthpop -h``.

The generated catalogs are saved at the location defined in the configuration (default: your_current_directory/output-files).


Import SynthPop to other script 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Importing SynthPop to another script allows for more flexibility. 
To do so, ensure that the parent directory is within your Python path and use the following code::
  
  import synthpop
  
  model = synthpop.SynthPop(config_file_name, **kwargs)
  model.init_populations()
  
All attributes of the configuration can also be specified by a keyword argument. 
It is then possible to run all specified locations via::
  
  model.process_all() 
  
or to run a specified location only::

  data, distributions = model.process_location(
        l_deg, b_deg, solid_angle, save_data=True) 
  
While ``process_all()`` only saves the results to disk, ``process_location()`` also returns the dataframe and several distributions (currently only distance distributions for each population). The ``save_data`` flag can be set to false to only return the data without saving it as a file.
