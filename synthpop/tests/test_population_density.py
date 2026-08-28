import synthpop.modules.population_density.constant as constant_density_module
from synthpop.modules.population_density import density_from_grid
import pdb
import numpy as np
import pytest

def check_precision(v_true, v_calc, prec):
	assert np.abs(v_calc-v_true)/v_true < prec, \
			f"Required precision {prec} not met for "+ \
			f"value {v_true} estimate of {v_calc}."
	print(v_true, v_calc)

def test_mass_integration_constant():
	"""
	Test the integration used to determine total stellar mass of a population in a given
	field. Test two different field sizes each for each shape, using a constant density.
	"""
	constant_density = constant_density_module.Constant()
	# Small circular field
	field_scale_deg = 0.01
	constant_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='circle', 
						field_scale_deg=field_scale_deg, max_distance=25)
	total_mass_integ = constant_density.total_mass
	total_mass_ana = constant_density.rho * 4./3*np.pi*constant_density.max_distance**3 * \
		  np.pi*field_scale_deg**2*(np.pi/180)**2/(4*np.pi)
	check_precision(total_mass_ana, total_mass_integ, 1e-4)

	# Larger circular field
	field_scale_deg = 1.0
	constant_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='circle', 
						field_scale_deg=field_scale_deg, max_distance=25)
	total_mass_integ = constant_density.total_mass
	total_mass_ana = constant_density.rho * 4./3*np.pi*constant_density.max_distance**3 * \
		  np.pi*field_scale_deg**2*(np.pi/180)**2/(4*np.pi)
	check_precision(total_mass_ana, total_mass_integ, 1e-4)

	# Small square field
	field_scale_deg = 0.02
	constant_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='box', 
						field_scale_deg=field_scale_deg, max_distance=25)
	total_mass_integ = constant_density.total_mass
	total_mass_ana = constant_density.rho * 4./3*np.pi*constant_density.max_distance**3 * \
		  field_scale_deg**2*(np.pi/180)**2/(4*np.pi)
	check_precision(total_mass_ana, total_mass_integ, 1e-4)

	# Larger square field
	field_scale_deg = 2.0
	constant_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='box', 
						field_scale_deg=field_scale_deg, max_distance=25)
	total_mass_integ = constant_density.total_mass
	total_mass_ana = constant_density.rho * 4./3*np.pi*constant_density.max_distance**3 * \
		  field_scale_deg**2*(np.pi/180)**2/(4*np.pi)
	check_precision(total_mass_ana, total_mass_integ, 1e-4)

	# Small rectangular field
	field_scale_deg = [0.02,0.01]
	constant_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='box', 
						field_scale_deg=field_scale_deg, max_distance=25)
	total_mass_integ = constant_density.total_mass
	total_mass_ana = constant_density.rho * 4./3*np.pi*constant_density.max_distance**3 * \
		  field_scale_deg[0]*field_scale_deg[1]*(np.pi/180)**2/(4*np.pi)
	check_precision(total_mass_ana, total_mass_integ, 1e-4)

	# Larger rectangular field
	field_scale_deg = [2.0,1.0]
	constant_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='box', 
						field_scale_deg=field_scale_deg, max_distance=25)
	total_mass_integ = constant_density.total_mass
	total_mass_ana = constant_density.rho * 4./3*np.pi*constant_density.max_distance**3 * \
		  field_scale_deg[0]*field_scale_deg[1]*(np.pi/180)**2/(4*np.pi)
	check_precision(total_mass_ana, total_mass_integ, 1e-4)

def test_mass_grid_nsc():
    """
    Test mass integration for non-constant density where high resolution is required.
    We use a NSC as it changes on small scales and is easily encapsulated in full in 1 field.
    """
    nsc_density = density_from_grid.DensityFromGrid(moment_file='vasiliev2026_huston_nsc_grid.dat')
    total_mass_from_agama = 14430946.610300452
    
    # Try all 3 field shapes and required 0.1% precision
    field_scale_deg = 0.2
    nsc_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='circle', 
                                field_scale_deg=field_scale_deg, max_distance=25)
    total_mass_integ = nsc_density.total_mass
    check_precision(total_mass_from_agama, total_mass_integ, 1e-3)
    
    field_scale_deg = 0.3
    nsc_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='box', 
                                field_scale_deg=field_scale_deg, max_distance=25)
    total_mass_integ = nsc_density.total_mass
    check_precision(total_mass_from_agama, total_mass_integ, 1e-3)
    
    field_scale_deg = [0.3, 0.4]
    nsc_density.update_location(l_deg=0.0, b_deg=0.0, field_shape='box', 
                                field_scale_deg=field_scale_deg, max_distance=25)
    total_mass_integ = nsc_density.total_mass
    check_precision(total_mass_from_agama, total_mass_integ, 1e-3)
