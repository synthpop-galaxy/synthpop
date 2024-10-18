from ._extinction import ExtinctionLaw

class Fitzpatrick2009(ExtinctionLaw):
    def __init__(self, **kwargs, ):
        self.extinction_law_name = 'Fitzpatrick2009'
        self.law_ref_wavelength = 0.549

    def Alambda_Aref(self, eff_wavelength: float, R_V: float = 3.0,
                        alpha: float = 2.5) -> float:
            """
            Given an effective wavelength lambda_eff, calculate the relative extinction A_lambda/A_V

            Parameters
            ----------
            eff_wavelength : float
                Effective Wavelength of the filter for which the extinction should be determined.
                in micrometer
            R_V : float
                interstellar reddening parameter
            """

            k = (0.349 + 2.087*R_V) * (1.0 / (1.0 + (eff_wavelength / 0.507)**alpha)) - R_V

            #return Al_AV
            return (k / R_V) + 1.

