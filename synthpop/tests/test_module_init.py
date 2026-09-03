import synthpop
from glob import glob
import numpy as np
import pytest
#from synthpop.modules.age import Age
#from synthpop.populations

def test_age():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/age/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.age.Age,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        print(class_name, 'OK')

def test_metallicity():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/metallicity/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.metallicity.Metallicity,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
        print(class_name, 'OK')

def test_imf():
    mod_files = glob(synthpop.constants.SYNTHPOP_DIR +
                          '/modules/initial_mass_function/*.py')
    for f in mod_files:
        class_name = f.split('/')[-1][:-3]
        if class_name.startswith('_'): continue
        cl = synthpop.synthpop_utils.get_subclass(synthpop.modules.initial_mass_function.InitialMassFunction,
                    synthpop.synthpop_utils.ModuleKwargs.parse_obj({"name":class_name}),
                    initialize=True)
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
        cl.density(np.array([1.0]), np.array([0.0]), np.array([0.0]))
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
        print(class_name, 'OK')
