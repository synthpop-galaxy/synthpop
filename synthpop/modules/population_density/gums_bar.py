"""
Bar density model from the Gaia Universe Model Snapshot
"""

__all__ = ["GumsBar", ]
__date__ = "2023-04-03"
__author__ = "J. Klüter"

import numpy as np
import scipy.special
from ._population_density import PopulationDensity

class GumsBar(PopulationDensity):
    """
    Bar density model from the Gaia Universe Model Snapshot
    
    Attributes
    ----------
    n0 : float
    x0 : float
    y0 : float
    z0 : float
    sigma_cut_of : float
    c_para : float
    c_perp : float
    dz_bone : float
    x_bone : float
    r_max : float
    """

    def __init__(self, n0, x0, y0, z0, alpha, beta, gamma, c_perp, c_para,
            dz_bone=0, x_bone=0, r_max=np.inf, sigma_cut_of=1e-10, **kwargs):
        super().__init__(**kwargs)
        self.density_unit = 'number'
        self.n0 = n0
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.sigma_cut_of = sigma_cut_of

        self.c_para = c_para
        self.c_perp = c_perp
        self.dz_bone = dz_bone
        self.x_bone = x_bone
        self.r_max = r_max

        sa, sb, sg = np.sin([alpha, beta, gamma])
        ca, cb, cg = np.cos([alpha, beta, gamma])
        rx_gamma = np.array([[1, 0, 0], [0, cg, -sg], [0, sg, cg]]).T
        ry_beta = np.array([[cb, 0, -sb], [0, 1, 0], [sb, 0, cb]]).T
        rz_alpha = np.array([[ca, sa, 0], [-sa, ca, 0], [0, 0, 1]]).T

        self.rot = np.matmul(np.matmul(rx_gamma, ry_beta), rz_alpha)

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
        n : ndarray [#/kpc^3]
            density at the given location in number density.

        """
        xyz_kpc = self.coord_trans.rphiz_to_xyz(r, phi_rad, z)
        x, y, z = np.dot(self.rot, xyz_kpc)

        r_perp = (np.abs(x / self.x0) ** self.c_perp
                  + np.abs(y / self.y0) ** self.c_perp) ** (1 / self.c_perp)

        # the documentation gives z_po = self.z0 * (1 + self.dz_bone * np.sin(x / self.x_bone)),
        # However when using abs(x), the density distribution aligns with the data of the
        # DR3 Gaia Universe Model
        z_po = self.z0 * (1 + self.dz_bone * np.sin(abs(x) / self.x_bone))

        rs = (r_perp ** self.c_para + np.abs(z / z_po) ** self.c_para) ** (1 / self.c_para)
        n = self.n0 * np.cosh(rs)**(-2)
        n *= np.exp(-0.5*np.maximum(np.sqrt(x**2+y**2)-self.r_max, 0)**2/self.sigma_cut_of**2)

        return n
