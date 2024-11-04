import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

import numpy as np

sys.path.append('..')
from synthpop_main import SynthPop

x = np.linspace(15, -15, 1000)
y = np.linspace(15, -15, 1000)
Y, X = np.meshgrid(y, x)


def plot_density(data, title):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(title)
    im = ax.imshow(data, extent=(y[0], y[-1], x[-1], x[0]))
    cbar = fig.colorbar(im, ax=ax)
    ax.set_xlabel('Y')
    ax.set_ylabel('X')



def main(config):
    plot_density(X, "validata X-axis\n increases upwards")
    plot_density(Y, "validata Y-axis\n increases to the left")
    mod = SynthPop(config)
    mod.init_populations()

    for population in mod.populations:
        r,phi,z = population.coord_trans.xyz_to_rphiz(X,Y,2)
        r = r.reshape(-1)
        phi = phi.reshape(-1)
        z = z.reshape(-1)
        dens = population.population_density.density(r, phi, z).reshape(X.shape)
        plot_density(dens, population.name)

    plt.show()

if __name__ == '__main__':
    if len(sys.argv) == 1:
        conf = 'gaia_universe_model.synthpop_conf'
    else:
        conf = sys.argv[1]
    main(conf)