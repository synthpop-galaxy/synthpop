__all__ = ["SpiseaImf", ]
__date__ = "2024-11-20"

import numpy as np
from spisea.imf import imf, multiplicity

try:
    from ._initial_mass_function import InitialMassFunction
except ImportError:
    from _initial_mass_function import InitialMassFunction

from typing import Callable

class SpiseaImf(InitialMassFunction):
    """
    Import an IMF module from SPISEA
    spisea_imf_name options are: IMF_broken_powerlaw, IMFSalpeter1955, Miller_Scalo_1979, Kennicutt_1983, 
        Kroupa_2001, Weidner_Kroupa_2004, 
    """

    def __init__(
            self, spisea_imf_name="Kroupa_2001",
            min_mass=0.1, max_mass=150, add_companions=False, multiplicity_kwargs=None, **kwargs
            ):
        """
        Initialize the IMF class for a Population class
        """
        spisea_imf_class=getattr(imf,spisea_imf_name)
        mass_limits = np.array([min_mass,max_mass])
        if not add_companions:
            self.spisea_imf = spisea_imf_class(massLimits=mass_limits)
        elif multiplicity_kwargs is None:
            self.spisea_imf = spisea_imf_class(massLimits=mass_limits, mutliplicity=multiplicity.MultiplicityResolvedDK())
        else:
            self.spisea_imf = spisea_imf_class(massLimits=mass_limits, mutliplicity=multiplicity.MultiplicityResolvedDK(**multiplicity_kwargs))
        self.add_companions = add_companions

    def imf(self, m_in):
        """
        empty function - not for use
        """
        return