# Script to run time tests for SynthPop vs. Galaxia
# Time test results saved to 'sp_gal_time_tests.txt'

import subprocess
import numpy as np
import synthpop
import pandas as pd
import time
from popsycle.synthetic import write_galaxia_params

# Test direction, arbitrary bulge direction
lon,lat = -1.4,1
sol_angs = [1e-5,1e-4,1e-3,1e-2,1e-1]

# SynthPop initialization
t0 = time.time()
model = synthpop.SynthPop(model_name='besancon_Robin2003',name_for_output='Besancon_validation',obsmag=False, 
                        maglim=['2MASS_Ks', 999999, 'keep'], extinction_map_kwargs={'name':'galaxia_3d'}, 
                        solid_angle_unit='deg^2',
                        max_distance=25, overwrite=True, output_file_type='hdf5')
model.init_populations()
t1 = time.time()
sp_init_time = t1-t0

# SynthPop run loop
times = []
stars = []
times.append(time.time())
for sa in sol_angs:
    cat, _ = model.process_location(lon,lat, sa)
    times.append(time.time())
    stars.append(len(cat))

# Galaxia run loop
output_root = 'gtt'
galaxia_galaxy_model_filename = '/g/lu/code/galaxia/docs/galaxyModelParams_PopSyCLEv3.txt'
gtimes = []
gstars = []
for sa in sol_angs:
    write_galaxia_params(output_root=output_root,
                         longitude=lon,
                         latitude=lat,
                         area=sa,
                         output_dir='galaxia_time_tests',
                         r_max=25)
    
    # Execute Galaxia
    cmd = 'galaxia -r %s_galaxia_params.txt %s' % (output_root, galaxia_galaxy_model_filename)
    gtimes.append(time.time())
    gresult = subprocess.run(cmd.split(' '),capture_output=True)
    gtimes.append(time.time())
    gstars.append(int([x for x in gresult.stdout.decode("utf-8").split('\n')[-10].split(' ') if x][-1]))
    

# Arrange and save runtimes
time_data = pd.DataFrame(data={'solid_angle': np.append(0.0,sol_angs),
                   'synthpop_time': np.append(0.0, np.diff(times))+sp_init_time,
                   'galaxia_time': np.append(0.0,np.diff(gtimes)[0::2]),
                   'synthpop_stars': np.append(0, stars),
                   'galaxia_stars': np.append(0, gstars)})
print(time_data)
time_data.to_csv('sp_gal_time_tests.txt', index=False)
