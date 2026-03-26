"""
Single value IMF (made for PBHs)
"""

__all__ = ['SingleValue', ]

import numpy as np
from scipy.special import erf, erfinv

try:
    from ._initial_mass_function import InitialMassFunction
except ImportError:
    from _initial_mass_function import InitialMassFunction


class SingleValue(InitialMassFunction):
    """
    Single value initial mass function (IMF)
    """

    def __init__(
            self, mass: float, **kwargs
            ):
        """
        Parameters
        ----------
        mass : float [Msun]
            mass value for all objects
        """
        super().__init__(min_mass, max_mass)
        self.imf_name = 'SingleValue'

        # setup control parameters
        self.mass = mass  

    def imf(self, m_in):
        """ initial mass function """
        if not isinstance(m_in, np.ndarray):
            m = np.array([m_in])
        else:
            m = m_in
        # 0.4342.. == 1/ln(10)
        prob = m==self.mass
        if not isinstance(m_in, np.ndarray):
            return prob[0]
        return prob

    def average_mass(self,
            min_mass: Union[np.ndarray, float, None] = None,
            max_mass: Union[float, None] = None
            ) -> float:
        return self.mass

    def draw_random_mass(
            self,
            min_mass: Union[np.ndarray, float, None] = None,
            max_mass: Union[float, None] = None,
            N: Union[float, None] = None
            ) -> Union[np.ndarray, float]:

        return np.ones(N)*self.mass  
