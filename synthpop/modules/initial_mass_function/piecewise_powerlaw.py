"""
Initial mass function for a piecewise power law, e.g.:
    
        For M ≤ m1:    ξ(M) ∝ M^-a0
        
        For m1 < M ≤ m2:    ξ(M) ∝ M^-a1
        
        For m2 < M:           ξ(M) ∝ M^-a2
        
The number of splitpoints is modifiable.
"""

__all__ = ["PiecewisePowerlaw", ]
__date__ = "2022-06-28"

import numpy as np

try:
    from ._initial_mass_function import InitialMassFunction
except ImportError:
    from _initial_mass_function import InitialMassFunction

from typing import Callable

class PiecewisePowerlaw(InitialMassFunction):
    """
    Initial mass function generator for a piecewise power law
    
    Attributes:
    -----------
    min_mass : float
        lower initial mass limit
    max_mass : float
        upper initial mass limit
    alphas : ndarray [float]
        power law indices for the mass chunks from lower mass to higher
    splitpoints : ndarray [float]
        mass values where pieces split;
        length must be length alphas minus 1
    """

    def __init__(
            self, min_mass=None, max_mass=None,
            alphas: tuple[float] = (1), splitpoints: tuple[float] = (), **kwargs
            ):
        super().__init__(min_mass, max_mass)
        self.imf_name = 'Piecewise Powerlaw'

        self.alphas = alphas
        self.splitpoints = splitpoints
        (self.imf_parts, self.F_imf_parts, self.F_imf_inverse_parts
        ) = self.get_functions(alphas, splitpoints)

    def get_functions(self, alphas, splitpoints):
        Fsplitpoints = []
        As = [1]
        F_parts = []
        F_inv_parts = []
        imf_parts = []
        mass = self.min_mass
        mass0 = self.min_mass

        F0 = 0

        for i, mass in enumerate(splitpoints):
            imf_parts.append(self.get_imf_parts(mass, mass0, As[-1], alphas[i]))

            func = self.get_F_parts(mass, mass0, As[-1], alphas[i])
            F_parts.append(func)
            Fsplitpoints.append(func(mass) + F0)

            F_inv_parts.append(self.get_F_inv_parts(mass, mass0,
                Fsplitpoints[-1], F0, As[-1], alphas[i]))

            As.append(As[-1] * mass ** (alphas[i + 1] - alphas[i]))
            F0 = Fsplitpoints[-1]
            mass0 = mass

        F_parts.append(self.get_F_parts(self.max_mass, mass, As[-1], alphas[-1]))
        imf_parts.append(self.get_imf_parts(self.max_mass, mass, As[-1], alphas[-1]))
        F_inv_parts.append(self.get_F_inv_parts(self.max_mass, mass, np.inf,
            F0, As[-1], alphas[-1]))

        return imf_parts, F_parts, F_inv_parts

    def get_imf_parts(self, m_upper, m_lower, a, alpha):
        m_lower = np.maximum(m_lower, self.min_mass)

        def func(m):
            return a * m ** (-alpha) * (m >= m_lower) * (m < m_upper)

        return func

    def get_F_parts(self, m_upper, m_lower, a, alpha) -> Callable:
        m_lower = np.maximum(m_lower, self.min_mass)
        if alpha == 1:
            def func(m):
                val = a * np.log(np.minimum(m, m_upper)) - a * np.log(m_lower)
                return np.maximum(val, 0)
        else:
            def func(m):
                val = a / (1 - alpha) * np.minimum(m, m_upper) ** (1 - alpha) \
                      - a / (1 - alpha) * m_lower ** (1 - alpha)
                return np.maximum(val, 0)
        return func

    def get_F_inv_parts(self, m_upper, m_lower, p_upper, p_lower, a, alpha):
        m_lower = np.maximum(m_lower, self.min_mass)
        if alpha == 1:
            def func(p):
                m = np.exp((p - p_lower) / a) * m_lower
                return m * (p > p_lower) * (p <= p_upper)
        else:
            def func(p):
                m = np.zeros(p.shape)
                gg = (p >= p_lower) & (p < p_upper)
                m[gg] = ((p[gg] - p_lower) * (1 - alpha) / a + m_lower ** (1 - alpha)) ** (
                            1 / (1 - alpha))
                return m
        return func

    # returns the number of stars of that initial mass
    def imf(self, m_in):
        """
        Initial mass function

        Parameters
        ----------
        m_in: initial mass

        Returns
        -------
        prob: probability at the initial mass

        """
        if not isinstance(m_in, np.ndarray):
            m = np.array([m_in])
        else:
            m = m_in

        prob = np.sum([i(m) for i in self.imf_parts], axis=0)

        if not isinstance(m_in, np.ndarray):
            return prob[0]
        return prob

    # The following part is optional, it speeds up the draw_random_mass
    # integration of imf

    def F_imf(self, m_in):
        """
        Integral of the initial mass function from 0 to m_in

        Parameters
        ----------
        m_in: initial mass

        Returns
        -------
        F: integral

        """
        if not isinstance(m_in, np.ndarray):
            m = np.array([m_in])
        else:
            m = m_in

        F = np.sum([f(m) for f in self.F_imf_parts], axis=0)

        if not isinstance(m_in, np.ndarray):
            return F[0]
        return F

    # inverse of F_imf
    def F_imf_inverse(self, p_in):
        """
        Inverse of F_imf

        Parameters
        ----------
        p_in : cumulative likelihood

        Returns
        -------
        met
        """

        if not isinstance(p_in, np.ndarray):
            p = np.array([p_in])
        else:
            p = p_in

        m = np.sum([f(p) for f in self.F_imf_inverse_parts], axis=0)

        if not isinstance(p_in, np.ndarray):
            return m[0]
        return m

    grid_max = 10000
    grid_min = 0
