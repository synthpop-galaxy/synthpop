# Quickstart for SynthPOP


## installation
  to install SynthPop, download the repository
  If you plan to import SynthPop to other scrips 
  ensure that the parent dictionary is in your Python Path.

### requirements
  The following packages are needed: 
  - numpy
  - scipy
  - pandas
  - request
  - json
  - tqdm
  - tables
  - dustmaps
  - astropy
  
  tables dustmaps and astropy are optional but highly recommended.
  install all by using
  ```pip install -r requirements.txt```

## run synthpop
  Synthpop can either be run as an individual script 
  or can be imported to other scripts

### run SynthPop as individual script
  to run the SynthPop type:
  ```
  python synthpop config_file 
  ```
  this process all locations as defined in the config_file 
  Have a look at ```python synthpop -h ``` for additional arguments
  all results are stored in the defined output location
  (default: synthpop/output-files)
  
### import SynthPop to other script 
  SynthPop can be imported to another script, which provides more flexibility. 
  to do so, ensure that the parent directory is within your python Path.
  ```
  import synthpop
  
  model = synthpop.SynthPop('config_file',**kwargs)
  model.init_population()
  ```
  all attributes in the configuration can also modify by a keyword argument
  
  to run all specified locations
  ```
  model.process_all() 
  ```

  to run a specified locations
  ```
  data, distributions = model.process_location(
        l_deg, b_deg,solid_angle) 
  ```
  while ```process_all()``` only saves the results to disk,
  ```process_location()``` also returns the dataframe and several distributions
  (currently only distance distributions for each population)

## Configure SynthPop
  SynthPop is controlled by a config json file.
  The Default can be found in synthpop/config_files.
  To modify arguments, these must be saved in an additional
  config file and passed to the
  When an arguments are not specified, synthpop falls back to 
  the value in the default_config file. 

  Note that _default_config.json is sorted into different categories, 
  These can be (but don't need to be) translated into the config file. 
  Also note that arguments starting with an '#' are ignored. 
 

### select magnitudes 
  The ```"chosen_band"``` argument is used to specify the returned filters, 
  With the default MIST isochrone system, any available filter can be selected. 
  However, these have to follow the MIST column names (e.g. The PanSTARRS g magnitude can only specify by "PS_g")
  Note that the implement MIST module downloads the data automatically. 

