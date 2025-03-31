"""
Initial mass function from Chabrier (2003):
    
    For m< m_switch:
    
        xi = A/m/ln(10) * exp(- (lg(m) - lg(m_center))**2 / (2*sigma**2))
        
    For m>=m_switch
    
        xi = A * m**(-alpha)
"""

__all__ = ['Chabrier', ]

import numpy as np
from scipy.special import erf, erfinv

try:
    from ._initial_mass_function import InitialMassFunction
except ImportError:
    from _initial_mass_function import InitialMassFunction


class Chabrier(InitialMassFunction):
    """
    The initial mass function (IMF) from Chabrier 2003
    """

    def __init__(
            self, min_mass: float = None, max_mass: float = None, m_center=0.079, A_less=0.158,
            sigma=0.69, m_switch=1, alpha=2.3, **kwargs
            ):
        """
        Parameters
        ----------
        min_mass : float [Msun]
            min mass to be generated
        max_mass : float [Msun]
            max mass to be generated
        m_center:
            center_mass of the gauss distribution
        A_less:
            normalisation constant
            note that A_higher is calculated such that the imf is steady at m_switch.

        sigma: [log Msun]

        m_switch: float [Msun]
            mass where to switch to a power law
        alpha: float
            slope of the power law
        """
        super().__init__(min_mass, max_mass)
        self.imf_name = 'Chabrier'

        # setup control parameters
        self.A_less = A_less
        self.center = np.log10(m_center)
        self.sigma = sigma
        self.m_switch = m_switch
        self.x = alpha - 1

        # Ensure that the IMF is smooth at the connection
        self.A_higher = (self.A_less
                        * np.exp(-(np.log10(self.m_switch) - self.center) ** 2 / (2 * sigma ** 2))
                        * self.m_switch ** self.x)

        # Estimate the integral from 0 to m_switch
        # needed for F_imf_inverse works properly ..
        self.F_m_switch = self.A_less * self.gaussian_integral(m_switch)
    def gaussian_integral(self,m_in):
        integral = np.sqrt(2 * np.pi) * self.sigma / 2 * (
                erf((np.log10(m_in) - self.center) / (np.sqrt(2) * self.sigma)) + 1)
        return integral

    def inverted_gaussian_integral(self, prop):
        arg = 2 * prop / (np.sqrt(2 * np.pi) * self.sigma) - 1
        m = 10 ** (self.center + np.sqrt(2) * self.sigma * erfinv(arg))
        return m

    def imf(self, m_in):
        """ initial mass function """

        if not isinstance(m_in, np.ndarray):
            m = np.array([m_in])
        else:
            m = m_in

        # 0.4342.. == 1/ln(10)
        prob = np.zeros(m.shape)
        prob[m >= self.m_switch] = self.A_higher * 0.4342944819 \
                                   * m[m >= self.m_switch] ** (-self.x-1)
        prob[m < self.m_switch] = self.A_less * 0.4342944819 / m[m < self.m_switch] \
            * np.exp(-(np.log10(m[m < self.m_switch]) - self.center) ** 2 / (2 * self.sigma ** 2))

        if not isinstance(m_in, np.ndarray):
            return prob[0]
        return prob

    # The following part is optional, it speeds up the draw_random_mass
    # integration of imf

    def F_imf(self, m_in):
        """ integral of the imf """

        # translate floats into ndarrays
        if not isinstance(m_in, np.ndarray):
            m = np.array([m_in])
        else:
            m = m_in

        F = np.zeros(m.shape)
        F[m > self.m_switch] = self.A_higher * 0.4342944819/(-self.x) \
                               * (m[m > self.m_switch]**(-self.x) - self.m_switch**(-self.x)) \
                               + self.F_m_switch
        F[m < self.m_switch] = self.A_less * self.gaussian_integral(m[m < self.m_switch])


        if not isinstance(m_in, np.ndarray):
            # translate ndarrays into floats if necessary.
            return F[0]
        return F

    def F_imf_inverse(self, p_in):
        """ inverse of F_imf """
        # translate floats into ndarrays
        if not isinstance(p_in, np.ndarray):
            p = np.array([p_in])
        else:
            p = p_in

        m = np.zeros(p.shape)
        m[p >= self.F_m_switch] = (-self.x  /self.A_higher / 0.4342944819
                                * (p[p >= self.F_m_switch]-self.F_m_switch)
                                + self.m_switch**(-self.x))**(-1/self.x)

        m[p < self.F_m_switch] = self.inverted_gaussian_integral(
                p[p < self.F_m_switch] / self.A_less)

        if not isinstance(p_in, np.ndarray):
            # translate ndarrays into floats if necessary.
            return m[0]
        return m

    grid_max = 10000
    grid_min = 0.001
