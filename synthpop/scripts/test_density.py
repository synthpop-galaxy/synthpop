import sys
import os

sys.path.append(f'{os.path.dirname(__file__)}/..')
import numpy as np
import synthpop
import matplotlib.pyplot as plt
import pandas as ps

ct = synthpop.synthpop_utils.coordinates_transformation
# select random location 
solid_angle = 1e-4
loc = np.random.uniform([-90, 10], [90, 45])
if loc[0] < 0: loc[0] += 360
# Generate field 
mod = synthpop.SynthPop(obsmag=False, maglim=['2MASS_Ks', 999999, 'keep'])
mod.init_populations()

w = 0.5
dd = np.arange(0, mod.parms.max_distance + w, w)
ddc = (dd[:-1] + dd[1:]) / 2
r_phi_z = ct.dlb_to_rphiz(dd, *loc)
pop_index = [i for i, pp in enumerate(mod.populations) if
    pp.population_density.density_unit != 'number']
dens = np.array([pp.population_density.density(*r_phi_z) for pp in mod.populations])
volume = (w * ddc ** 2 + w ** 3 / 12) * solid_angle

data, cc = mod.process_location(*loc, solid_angle, save_data=False)
bins = ps.cut(data.Dist, dd)
bined_data = data.groupby(bins)
color = iter(plt.cm.rainbow(np.linspace(0, 1, len(pop_index))))

for i in pop_index:
    c = next(color)
    plt.semilogy(dd, dens[i], ':', c=c)
    binned_data_i = data.where(data['pop'] == i).groupby(bins)
    masses_i = binned_data_i.Mass.sum()
    plt.semilogy(ddc[masses_i > 1], masses_i[masses_i > 1] / volume[masses_i > 1], 'x-', c=c)

dens_tot = np.sum(dens, axis=0)
plt.semilogy(dd, dens_tot, 'grey', lw=3)

binned_data = data.groupby(bins)
masses = binned_data.Mass.sum()
plt.semilogy(ddc, masses / volume, 'kx')
plt.ylim([0.5 / volume[-1], 10 * max(masses / volume)])
plt.xlabel("dist [kpc]")
plt.ylabel("density [M_sun/kpc^3]")
plt.show()
