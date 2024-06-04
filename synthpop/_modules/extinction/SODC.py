from ._extinction import ExtinctionLaw

alpha = 2.255

class SODC(ExtinctionLaw):
	def __init__(self, **kwargs, ):
		self.extinction_law_name = 'NewExtLaw'

	def Alambda_AV(self, eff_wavelength: float, R_V: float = 3.1) -> float:
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

        	x = 1. / eff_wavelength
        	if x >= 1.1:
        		y = x - 1.82
        		a = 1 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4 - 1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8
        		b = 1.952*y + 2.908*y**2 - 3.989*y**3 - 7.985*y**4 + 11.102*y**5 + 5.491*y**6 - 10.805*y**7 + 3.347*y**8
        		#Al_AV = a+b/R_V
        	

        	else:
        		xpow = x ** 2.255
        		a = 0.564 * xpow
        		b = -0.484 * xpow
        		#Al_AV = a*x**alpha+b*x**alpha/R_V
        	
        	#return Al_AV
        	return a + b / R_V
