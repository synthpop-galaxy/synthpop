"""
This file includes Sun class
which stores the location and velocity of the sun and the velocities of the lsr
"""

__all__ = ["SunInfo", "default_sun"]
__author__ = "J. Klüter, S. Johnson, M.J. Huston"
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__data__ = "2023-05-05"

from functools import cached_property
import numpy as np
from pydantic import BaseModel

class SunInfo(BaseModel):
    # location of the sun
    x: float = -8.178
    y: float = 0.0
    z: float = 0.017

    # rectangular motion of the sun
    # from Reid & Brunthaler (2020)
    u: float = 12.9  # [km/s] = -V_R
    v: float = 245.6  # [km/s] = V_phi
    w: float = 7.78  # [km/s] = V_z

    # velocity of the local standard of rest
    # from Schönrich et al. (2010)
    u_lsr: float = 1.8
    v_lsr: float = 233.4
    w_lsr: float = 0.53

    # degree direction of the solar apex in galactic coordinates
    l_apex_deg: float = 53.
    b_apex_deg: float = 25.

    class Config():
        try: #pydantic version compatibility
            keep_untouched = (cached_property,)
        except:
            ignored_types = (cached_property,)

    @cached_property
    def gal_dist(self) -> float:
        # distance to the galactic center (kpc)
        return (self.x ** 2 + self.y ** 2 + self.z ** 2) ** 0.5

    @cached_property
    def r(self) -> float:
        # distance to the galactic center (kpc) projected on the galactic plane
        return (self.x ** 2 + self.y ** 2) ** 0.5

    @cached_property
    def theta(self) -> float:
        return np.arcsin(self.z / self.gal_dist)

    @cached_property
    def cos_theta(self) -> float:
        return np.cos(self.theta)

    @cached_property
    def sin_theta(self) -> float:
        return np.sin(self.theta)

    @cached_property
    def v_pec(self) -> float:
        # peculiar motion of the Sun relative to the LSR
        return ((self.u - self.u_lsr) ** 2
                 + (self.v - self.v_lsr) ** 2
                 + (self.w - self.w_lsr) ** 2)**0.5

default_sun = SunInfo()
