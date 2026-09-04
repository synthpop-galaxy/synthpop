import synthpop

# Make sure we can initialize all the existing models
# We swap in a non-default extinction map for loading speed

def test_init_besancon_Robin2003():
    mod = synthpop.SynthPop(model_name='besancon_Robin2003',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_galaxia_Sharma2011():
    mod = synthpop.SynthPop(model_name='galaxia_Sharma2011',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_genstars_Koshimoto2021():
    mod = synthpop.SynthPop(model_name='genstars_Koshimoto2021',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_genstars_Koshimoto2022():
    mod = synthpop.SynthPop(model_name='genstars_Koshimoto2022',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_GUMS_dr3():
    mod = synthpop.SynthPop(model_name='GUMS_dr3',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_GUMS_dr3_mod_dens():
    mod = synthpop.SynthPop(model_name='GUMS_dr3_mod_dens',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_Huston2024():
    mod = synthpop.SynthPop(model_name='Huston2024',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_Huston2025():
    mod = synthpop.SynthPop(model_name='Huston2025',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_Huston2026():
    mod = synthpop.SynthPop(model_name='Huston2026',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()

def test_init_validation_model():
    mod = synthpop.SynthPop(model_name='validation_model',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()
    
def test_init_Vasiliev2026():
    mod = synthpop.SynthPop(model_name='Vasiliev2026',
                         extinction_map_kwargs={"name":"maps_from_dustmaps", "dustmap_name":"marshall"},
                        )
    mod.init_populations()