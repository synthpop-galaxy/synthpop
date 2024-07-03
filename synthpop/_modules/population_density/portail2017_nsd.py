"""
Nuclear stellar disk
"""

__date__ = "2021-10-18"
__all__ = ["nsd", ]

import numpy as np
import scipy.special
from ._population_density import PopulationDensity


class portail(PopulationDensity):
    def __init__(
            self, rho0=2.0e9*51, h_r = 0.250, h_z = 0.050, bar_angle=29.4,
            **kwargs
            ):
        """
        NSD model in testing
        """
        super().__init__(**kwargs)
        self.rho0 = rho0
        self.h_r = h_r
        self.h_z = h_z
        self.bar_angle = bar_angle*np.pi/180
        self.density_unit = 'mass'

    def density(self, r, phi_rad, z):
        alpha = np.pi/2.0 - self.bar_angle

        # Rotate into bar frame
        xb = -r * np.cos(phi_rad - self.bar_angle)
        yb = r * np.sin(phi_rad - self.bar_angle)

        return self.rho0 * np.exp(-np.sqrt(xb**2 + (yb/0.5)**2)/self.h_r) * np.exp(-np.abs(z)/self.h_z)
