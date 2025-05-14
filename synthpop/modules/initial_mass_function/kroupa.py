"""
Initial mass function from Kroupa (2001):
    
    For 0.01 < M ≤ 0.08:    ξ(M) ∝ M^(-0.3)
    
    For 0.08 < M ≤ 0.50:    ξ(M) ∝ M^(-1.3)
    
    For 0.50 < M:           ξ(M) ∝ M^(-2.3)
"""

__all__ = ["Kroupa", ]
__date__ = "2022-06-28"

import numpy as np
from ._initial_mass_function import InitialMassFunction

class Kroupa(InitialMassFunction):
    """
    IMF from Kroupa (2001)
    """

    def __init__(self, min_mass=None, max_mass=None, **kwargs):
        """
        Initialize the IMF class for a Population class
        """
        super().__init__(min_mass, max_mass)
        self.imf_name = 'Kroupa'


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

        prob = 1.78513 * m ** -0.3
        prob[m > 0.08] = 0.142810 * m[m > 0.08] ** -1.3
        prob[m > 0.5] = 0.0714053 * m[m > 0.5] ** -2.3
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
        F = 1.78513 / 0.7 * m ** 0.7
        i = m > 0.08
        F_008 = 0.4352460193816074  # integral from 0 to 0.08
        F[i] = -0.142810 / 0.3 * (m[i] ** -0.3 - 0.08 ** -0.3) + F_008
        i = m > 0.5
        F_05 = 0.8647514214666863  # integral from 0 to 0.5
        F[i] = -0.0714053 / 1.3 * (m[i] ** -1.3 - 0.5 ** -1.3) + F_05
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
        F_05 = 0.8647514214666863  # integral from 0 to 0.08
        F_008 = 0.4352460193816074  # integral from 0 to 0.5

        m = (p * 0.7 / 1.78513) ** (1 / 0.7)
        m[p > F_008] = (-(p[p > F_008] - F_008) * 0.3 / 0.142810 + 0.08 ** -0.3) ** (-1 / 0.3)
        m[p > F_05] = (-(p[p > F_05] - F_05) * 1.3 / 0.0714053 + 0.5 ** -1.3) ** (-1 / 1.3)
        if not isinstance(p_in, np.ndarray):
            return m[0]
        return m

    grid_max = 10000
    grid_min = 0
