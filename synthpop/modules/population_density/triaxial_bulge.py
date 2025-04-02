"""
Bulge density models for all triaxial formulations
from Dwek (1995).
"""

__all__ = ["TriaxialBulge", ]
__author__ = "M.J. Huston"
__date__ = "2024-02-15"

import numpy as np
import scipy.special
from .. import const
from ._population_density import PopulationDensity

class TriaxialBulge(PopulationDensity):
    """
    Bulge density models for all triaxial formulations from Dwek (1995).
    
    Attributes
    ----------
    triaxial_type : str
        functional form as defined in Dwek (1995);
        available options are 'G1', 'G2', 'G3', 'E1', 'E2', 'E3', 'P1', 'P2', 'P3'
    density_unit: str
        options: 'mass' or 'number'
    x0, y0, z0 : float [kpc]
        semi major axes of the bulge
    rho0 : float [m_sun/kpc^3 or stars/kpc^3]
        central density, unit determined by density_unit
    Rmax : float [kpc]
        cutoff radius for bulge
    bar_angle : float [degrees]
        angle of the bar
    """

    def __init__(self, triaxial_type: str, density_unit: str, x0: float, y0: float, z0: float, rho0: float, Rmax=np.inf, bar_angle=29.4, **kwargs):
        super().__init__()
        self.population_density_name = "Bulge_Density"
        self.density_unit = density_unit
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.rho0 = rho0
        self.Rmax = Rmax
        self.bar_ang = bar_angle * np.pi / 180

       #Functional form
        if triaxial_type not in ['G1','G2','G3','E1','E2','E3','P1','P2','P3']:
            raise ValueError("Invalid triaxial bulge type: "+triaxial_type)
        # Eqns 3a-i, 4a-c
        def rho_G1(x,y,z):
            rr = np.sqrt((x/self.x0)**2 + (y/self.y0)**2 + (z/self.z0)**2)
            return self.rho0 * np.exp(-0.5*rr**2)
        def rho_G2(x,y,z):
            rr_e = np.absolute(x)/self.x0 + np.absolute(y)/self.y0 + np.absolute(z)/self.z0
            return self.rho0 * np.exp(-0.5*rr_e**2)
        def rho_G3(x,y,z):
            rr = np.sqrt((x/self.x0)**2 + (y/self.y0)**2 + (z/self.z0)**2)
            return self.rho0 * rr**-1.8 * np.exp(-rr**3)
        def rho_E1(x,y,z):
            rr_e = np.absolute(x)/self.x0 + np.absolute(y)/self.y0 + np.absolute(z)/self.z0
            return self.rho0 * np.exo(-rr_e)
        def rho_E2(x,y,z):
            rr = np.sqrt((x/self.x0)**2 + (y/self.y0)**2 + (z/self.z0)**2)
            return self.rho0 * np.exp(-rr)
        def rho_E3(x,y,z):
            rr_s = (((x/self.x0)**2 + (y/self.y0)**2)**2 + (z/self.z0)**4)**0.25
            return self.rho0 * scipy.special.kn(0, rr_s)
        def rho_P1(x,y,z):
            rr = np.sqrt((x/self.x0)**2 + (y/self.y0)**2 + (z/self.z0)**2)
            return self.rho0 * 1/(1+r)**4       
        def rho_P2(x,y,z):
            rr = np.sqrt((x/self.x0)**2 + (y/self.y0)**2 + (z/self.z0)**2)
            return self.rho0 * 1/(r*(r+1)**3)   
        def rho_P3(x,y,z):
            rr = np.sqrt((x/self.x0)**2 + (y/self.y0)**2 + (z/self.z0)**2)
            return self.rho0 * 1/(1+r**2)**2
        # Select proper function
        rho_func_dict = {'G1':rho_G1, 'G2':rho_G2, 'G3':rho_G3, 'E1':rho_E1, 'E2':rho_E2, 'E3':rho_E3, 'P1':rho_P1, 'P2':rho_P2, 'P3':rho_P3}
        self.rho_func = rho_func_dict[triaxial_type]


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
        # Align coordinates with the bar,
        xb = -r * np.cos(phi_rad - self.bar_ang)
        yb = r * np.sin(phi_rad - self.bar_ang)
        zb = z

        # Apply cutoff radius (Eqn 5)
        # 1 when r<Rmax, exponential drop-off otherwise
        r_0 = 0.5
        f_cut = (r<=self.Rmax).astype(int) + (r>self.Rmax).astype(int)*np.exp(-r**2/2/r_0**2)

        return self.rho_func(xb,yb,zb)*f_cut
