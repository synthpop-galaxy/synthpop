"""
This script is used for the validation of the different modules
It creates the several plots showing the input distribution and the resulting distribution.
"""



__author__ = 'Jonas Klüter'
__date__ = '2023-06-29'

import sys
import os
import numpy as np
import pandas
import matplotlib.pyplot as plt
import pandas as pd

DIRNAME = os.path.abspath(os.path.dirname(__file__))

# import synthpop
import synthpop

# set up SynthPop model
model = synthpop.SynthPop(os.path.join(DIRNAME,'validation_config.json'))
# initialize populations
model.init_populations()

# place to collect the data from different line of sights
data = {}
# loop over all location
for loc in model.get_iter_loc():
    # run synthpop for the given location and solid angle
    data_loc, distribution = model.process_location(
        *loc, model.parms.solid_angle, model.parms.solid_angle_unit,
        save_data=False)

    # store data in the dictionary
    data[loc] = data_loc
    break
combined_data = pandas.concat([d for d in data.values()])


def get_color_iter(n):
    return iter(['red', 'green', 'blue', 'gold'])
def validate_imf(comp_data, populations=(0, 1, 2)):
    """ plot the given imf and resulting distribution"""
    color = get_color_iter(len(populations))
    plt.figure(figsize=(4,4))
    x = np.logspace(np.log10(model.parms.mass_lims['min_mass']),
                    np.log10(model.parms.mass_lims['max_mass']), 1001)

    for i in populations:
        c=next(color)

        # get the imf for the current population
        imf = model.populations[i].imf
        # normalise the IMF to an area of 1
        norm = imf.F_imf(model.parms.mass_lims['max_mass'])\
               -imf.F_imf(model.parms.mass_lims['min_mass'])
        # plot input imf
        plt.loglog(x, imf.imf(x)/norm,c=c, alpha=0.5, lw=3)

        # select stars generated for the current population
        population_mask = combined_data['pop'] == i
        # plot histogramm
        plt.hist(comp_data.iMass[population_mask], bins=x[::25],
            histtype='step', color=c, density=True, label=imf.imf_name)
    plt.title('IMF')
    plt.ylabel('probability density')
    plt.xlabel('initial mass [M$_{\odot}]$')
    plt.legend()

def validate_age(comp_data, populations=(0, 1, 2)):
    """ plot the given imf and resulting distribution"""
    color = get_color_iter(len(populations))
    plt.figure(figsize=(4,4))
    x = np.linspace(0,14, 1001)

    for i in populations:
        c=next(color)

        # get the imf for the current population
        age = model.populations[i].age
        # normalise the IMF to an area of 1

        # plot input metallicity
        plt.hist(age.draw_random_age(10_000_000), bins=x, density=True,
            histtype="step", lw=3, color=c, alpha=0.5)

        #select stars generated for the current population
        population_mask = combined_data['pop'] == i
        #plot histogramm
        plt.hist(comp_data.loc[population_mask,'age'], bins=x[::25],
            histtype='step', color=c, density=True, label=age.age_func_name)
    plt.title('Age')
    plt.ylabel('probability density')
    plt.xlabel('age [Gyr]$')
    plt.legend()

def validate_met(comp_data, populations=(0, 1, 2)):
    """ plot the given imf and resulting distribution"""
    color = get_color_iter(len(populations))
    plt.figure(figsize=(4,4))
    x = np.linspace(-3,1, 1001)

    for i in populations:
        c=next(color)

        # get the imf for the current population
        met = model.populations[i].metallicity
        # normalise the IMF to an area of 1
        plt.plot(x, met.likelyhood_distribution(x),c=c,lw=3, alpha=0.5)

        #select stars generated for the current population
        population_mask = combined_data['pop'] == i
        #plot histogramm
        plt.hist(comp_data.loc[population_mask,'Fe/H_initial'], bins=x[::25],
            histtype='step', color=c, density=True, label=met.metallicity_func_name)

    plt.title('Metallicity')
    plt.ylabel('probability density')
    plt.xlabel('Fe/H [dex]')
    plt.legend()


validate_imf(combined_data)
validate_met(combined_data)
validate_age(combined_data)


plt.show()