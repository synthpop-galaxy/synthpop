#from position import Position
from initialmass import InitialMass
#from massdistributions import MassDist
from age import Age
from metallicity import Metallicity
#from kinematics import Kinematics
from evolution import Evolution
#from populationdensity import PopulationDensity
#from coordinates import Coordinates
from calibration_parameters import dMax,dStep,m_stars,barang,chosenBands, bandOpt,maglim,magsys,angleMethod,outputLocation,outputFileType,populations,l_set,b_set,x_sun,y_sun,z_sun,obsmag,columnWidth,cols,starsPerPop, solidAngle, mass_lims

from numpy import zeros, sqrt, arctan2, arccos, arcsin, sin, cos, tan, pi, log10,searchsorted, random, arange
#from astropy.utils.console import ProgressBar

import pandas
import time
import pickle
import pdb


def make_calib_file():
    for pop in populations:
        data = np.loadtxt('pop_'+str(pop)+'_nstar_100000.cal.csv',skiprows=1)
    return


def calibrate_pop_density(popid,n_stars=10000):
    """                                                                                                                                                                                                                                                                                                             
    A function that numerically determines the ratio of the initial mass to the final mass                                                                                                                                                                                                                          
    given a population                                                                                                                                                                                                                                                                                              
    """
    #set up isochrones
    metToFile = ['m4.00', 'm3.50', 'm3.00', 'm2.50', 'm2.00', 'm1.75', 'm1.50', 'm1.25', 'm1.00', 'm0.75', 'm0.50', 'm0.25', 'p0.00', 'p0.25', 'p0.50']
    #pull out all the information from the isochrone file and puts it into tracks under the index of the lowest
    isochrones = {}
    for met in metToFile:
        fname = 'mistIsochrones/' + magsys + '/MIST_v1.2_feh_' + met + '_afe_p0.0_vvcrit0.4_' + magsys + '.iso.cmd'
        #fname: name of the file
        #delim_whitespace=true: uses whitespaces as delimiting characters                                                                                                                                                                                                                                           
        #dtype=float: makes each data entry in the df type float                                                                                                                                                                                                                                                    
        #comment='#': ignores any text after an '#' on the line                                                                                                                                                                                                                                                     
        #skip_blank_lines=true: ignores any blank lines                                                                                                                                                                                                                                                             
        isochrones[met] = pandas.read_csv(fname,delim_whitespace=True,comment='#',skip_blank_lines=True,low_memory=False)
        
          
        

    #columns that will be reported
    headers = cols.copy()

    df = pandas.DataFrame(columns=headers)

    total_init_mass = 0
    total_final_mass = 0

    #need to be able to generate age, metalicity for the population
    im = InitialMass()
    a = Age( popid)
    metallicity = Metallicity(popid)
    evol = Evolution()
    bands = chosenBands

    start_time = time.process_time()
    for star in range(n_stars):

        m_initial = im.next_mass()
        s_age = a.rand_age()*10**9
        s_met = metallicity.rand_met()
        if m_initial > 0.1:
            initmag, s_props, idx = evol.magcut(isochrones, m_initial, s_met, s_age, bands)
            m_evolved = s_props['star_mass']
            log_l = s_props['log_L']
            log_teff = s_props['log_Teff']
            log_g = s_props['log_g']

        else:
            m_evolved = m_initial
            log_l = -10
            log_teff = -10
            log_g = -10
            
        total_init_mass += m_initial
        total_final_mass += m_evolved
        df.loc[star] = [m_initial,m_evolved,popid,s_met,log_l,log_teff,log_g,log10(s_age)]
        if star%1000 == 0:
            print('On star ',star,' of ',n_stars,' for population ',popid)
            print('Total init mass drawn: ',total_init_mass)
            print('Running init/final : ',total_init_mass/total_final_mass)
            print('Running final/init : ',total_final_mass/total_init_mass)

        #make sure to temp write every 10000
        if (star+1)%10001==0:
            #write output
            filename = outputLocation + '/pop_' + str(popid) + '_nstar_' + str(n_stars) + '.cal.' + outputFileType
            if outputFileType == 'pickle':
                df.to_pickle(filename)
            elif outputFileType == 'csv':
                df.to_csv (filename, index=False, header=True,float_format='%0.6e')
            else:
                print("ERROR: Invalid output file type")

            
    print(" time per star  " + str((time.process_time()-start_time)/n_stars))
    print("Took total time " + str(time.process_time()))
    print("Total initial mass: ",total_init_mass)
    print("Total final mass: ",total_final_mass)
    #write output
    filename = outputLocation + '/pop_' + str(popid) + '_nstar_' + str(n_stars) + '.cal.' + outputFileType
    if outputFileType == 'pickle':
        df.to_pickle(filename)
    elif outputFileType == 'csv':
        df.to_csv (filename, index=False, header=True, float_format='%0.6e')
    else:
        print("ERROR: Invalid output file type")


if __name__=='__main__':
    pdb.set_trace()
    for popid in populations:
        print('Population: ',popid)
        calibrate_pop_density(popid,n_stars=100000)
