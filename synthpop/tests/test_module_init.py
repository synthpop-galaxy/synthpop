import synthpop
from glob import glob
import numpy as np
import pytest
import pdb
#from synthpop.modules.age import Age
#from synthpop.populations

# Make sure we can initialize and draw from each Age module
def test_age():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/age/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.age.Age,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        age = cl.draw_random_age(1000)
        assert ~np.any(np.isnan(age))
        print(class_name, 'OK')

# Make sure we can initialize and draw from each Metallicity module
def test_metallicity():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/metallicity/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.metallicity.Metallicity,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        met = cl.draw_random_metallicity(1000, x=np.zeros(1000), y=np.zeros(1000), z=np.zeros(1000))
        assert ~np.any(np.isnan(met))
        print(class_name, 'OK')

# Make sure we can initialize and draw from each standard IMF module
def test_imf():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/initial_mass_function/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.initial_mass_function.InitialMassFunction,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        if class_name=='spisea_imf': continue
        m = cl.draw_random_mass(0.01,1e4,100)
        assert ~np.any(np.isnan(m))
        print(class_name, 'OK')

def test_multiplicity():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/multiplicity/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.multiplicity.Multiplicity,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        # TODO: test running it probably...
        print(class_name, 'OK')

def test_ifmr():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/initial_final_mass_relation/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.initial_final_mass_relation.InitialFinalMassRelation,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        # TODO: test running it probably...
        print(class_name, 'OK')

def test_density():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/population_density/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.population_density.PopulationDensity,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        dens = cl.density(np.linspace(0,10,100), np.zeros(100), np.zeros(100))
        assert ~np.any(np.isnan(dens))
        print(class_name, 'OK')

def test_kinematics():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/kinematics/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.kinematics.Kinematics,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        kine = cl.draw_random_velocity(np.linspace(0.001,10,100), np.zeros(100), np.zeros(100))
        assert ~np.any(np.isnan(kine))
        print(class_name, 'OK')

@pytest.mark.skip(reason="Slow")
def test_postproc():
    # make default model but with faster loading ext map
    model = synthpop.SynthPop(extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"})
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/post_processing/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        mod_kwargs = synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name})
        mod_kwargs.model = model
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.post_processing.PostProcessing,
                    mod_kwargs,
                    initialize=True)
        # TODO: run a test catalog through them
        print(class_name, 'OK')
