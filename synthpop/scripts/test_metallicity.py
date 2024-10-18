import matplotlib.pyplot as plt
import numpy as np

from . import double_gaussian as dg

ddd = {"weight": 0.323, "mean1": -0.31, "std1": 0.31, "mean2": 0.26, "std2": 0.20}

dd = dg.DoubleGaussian(**ddd)
ff = dd.draw_random_metallicity(10000)
print(ff)
plt.hist(ff, bins=10, density=True)
x = np.linspace(-1.5, 1, 100)
y1 = ddd['weight'] / ddd["std1"] / np.sqrt(2 * np.pi) * np.exp(
    -0.5 * ((x - ddd["mean1"]) / ddd["std1"]) ** 2)
y2 = (1 - ddd['weight']) / ddd["std2"] / np.sqrt(2 * np.pi) * np.exp(
    -0.5 * ((x - ddd["mean2"]) / ddd["std2"]) ** 2)
y = y1 + y2
plt.plot(x, y, 'k')
plt.plot(x, y1, 'r')
plt.plot(x, y2, 'r')

plt.show()
