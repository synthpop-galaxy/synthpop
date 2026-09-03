"""
Einasto density profiles, used e.g. for the thin disk in the
Besancon Model (Robin et al., 2003).
"""

__all__ = ["Besancon2003Thindisk", ]
__author__ = "J. Klüter"
__date__ = "2022-07-12"

import numpy as np
from .. import const
from ._population_density import PopulationDensity

class Besancon2003Thindisk(PopulationDensity):
    """
    Einasto profile

    Attributes
    ----------
    e : float
        disc axis ratio
    p0 : float [Msun/kpc**3]
        density at the position of sun
    hrp : float [kpc]
        scale radius of the positive Gaussian
    hrm : float [kpc]
        hrm is the scale radius of the disk hole
    offset : float
        modification of the exponent
    power :  float
        modification of the exponent
        for the Besancon thin disk:
        offset = 0; power = 2 for age<=0.15 Gyr
        offset = 0.5; power = 1 for age> 0.15 Gyr
    disk_cutoff : float [kpc]
        outer edge of the disc
    flare_flag : bool
        include flare
    """

    def __init__(
            self, e: float = 0.0140, p0: float = 0.003999e9,
            hrp: float = 5.0, hrm: float = 3.0,
            offset: float = 0, power: float = 1,
            disk_cutoff: float = 14, flare_flag: boolean = True,
            **kwargs
            ):
        super().__init__(**kwargs)
        self.population_density_name = "Besancon2003Thindisk"
        self.density_unit = 'mass'
        self.disk_cutoff = disk_cutoff
        self.p0 = p0
        self.e = e
        self.hrp = hrp
        self.hrm = hrm
        self.offset = offset
        self.power = power
        self.flare_flag = flare_flag

    def density(self, r: np.ndarray, phi_rad: np.ndarray, z: np.ndarray) -> np.ndarray:
        """

        Estimates the density at the given position

        Parameters
        ----------
        r : ndarray ['kpc']
            Distance to z axis
        phi_rad : ndarray ['rad']
            azimuth angle of the stars. phi_rad = 0 is pointing towards sun.
        z : height above the galactic plane (corrected for warp of the galaxy)

        Returns
        -------
        rho : ndarray [M_sun/kpc^3 or #/kpc^3]
            density at the given location, either in number density evolved
            mass density or initial mass density should be specified in density_unit.

        """
        if self.flare_flag:
            k_flare = self.get_kappa_flare(r)
            k_flare0 = self.get_kappa_flare(self.sun.r)
        else:
            k_flare = 1
            k_flare0 = 1

        a = np.sqrt(r ** 2 + (z / k_flare / self.e) ** 2)
        a0 = np.sqrt(self.sun.r ** 2 + (self.sun.z / k_flare0 / self.e) ** 2)

        def exp_arg(x):  # function to estimate the argument for exp()
            return -np.sqrt(self.offset ** 2 + x ** 2) ** self.power

        # density normalization at the Sun [solar mass per kpc^3]
        d0 = np.exp(exp_arg(a0 / self.hrp)) - np.exp(exp_arg(a0 / self.hrm))

        # estimate density
        rho = self.p0 / d0 * (np.exp(exp_arg(a / self.hrp)) - np.exp(exp_arg(a / self.hrm)))

        rho *= r <= self.disk_cutoff
        # radius of disk set to 0 if r > disk_cutoff. used this way so numpy arrays,
        # and floats works the same without checking the type

        return rho
