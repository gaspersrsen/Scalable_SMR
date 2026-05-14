from .geo import *

#TODO move outside initial builder

def init_mix(self):
    ### ENRICHMENT MIXING
    self.mix_fuel_mats = {}
    self.mix_fuel_mats["u235_O2"] = make_fuel_material(self, "U235_O2", 1)
    self.mix_fuel_mats["u238_O2"] = make_fuel_material(self, "U238_O2", 0)
    # print(u235_O2,u238_O2)
    the_fuel = self.fuel_materials[0]
    m_U = np.sum([the_fuel.get_mass_density(nuc) for nuc in ["U235", "U238"]])
    m_235 = np.sum([the_fuel.get_mass_density(nuc) for nuc in ["U235"]])
    # print(f"The fuel {m_235/m_U} enrichment")



def mat_mix_U(self, frac):#FINDER FUNCTION wo
    # init_mix(self)
    the_fuel = self.fuel_materials[0]
    u235_O2 = self.mix_fuel_mats["u235_O2"]
    u238_O2 = self.mix_fuel_mats["u238_O2"]
    mO16 = openmc.data.atomic_mass("O16")
    mU235 = openmc.data.atomic_mass("U235")
    mU238 = openmc.data.atomic_mass("U238")
    global mix_ao_U
    mats=[u235_O2, u238_O2]
    muu = 1/(frac/mU235 + (1-frac)/mU238)
    mix_ao_U, ao_fracs, wgts = openmc.search.get_ao_mix_materials(mats, [(frac*(muu/(2*mO16+muu))), None],
                                                                              fracs_target=["U235",None],
                                                                              percent_type='wo', return_wgts=True)#type can be anything
    the_fuel.update_material(mix_ao_U)
    m_mix = the_fuel.get_mass_density()
    m_U = np.sum([the_fuel.get_mass_density(nuc) for nuc in ["U235", "U238"]])
    m_235 = np.sum([the_fuel.get_mass_density(nuc) for nuc in ["U235"]])
    print(f"Rebuilt material with {m_235/m_U} enrichment and mass density {m_mix} g/cm3")
    out= {}
    out["materials"] = [the_fuel]
    print(np.array(list(ao_fracs.values())))
    out["nuc_fractions"] = [np.array(list(ao_fracs.values()))[:,0] - wgts[0]/wgts[1]*np.array(list(ao_fracs.values()))[:,1]]
    return out




### BORON MIXING
def mat_mix_B(self, conc):#FINDER FUNCTION wo
    global mix_ao
    pure_water = self.mats["pure_water"]
    pure_boric_acid = self.mats["pure_boric_acid"]
    mix_ao, ao_fracs, wgts = openmc.search.get_ao_mix_materials(
        [pure_water, pure_boric_acid], [None, conc*1e-6],
        fracs_target=[None, "B"],
        percent_type='wo', return_wgts=True)#type can be anything
    self.mats["cool"].update_material(mix_ao)
    m_mix = self.mats["cool"].get_mass_density()
    m_boron = np.sum([self.mats["cool"].get_mass_density(nuc) for nuc in ["B10","B11"]])
    print(f"Rebuilt material with {m_boron/m_mix*1e6} ppm")
    out= {}
    out["materials"] = [self.mats["cool"]]
    # out["nuc_fractions"] = [np.array(list(ao_fracs.values()))[:,1]]
    out["nuc_fractions"] = [np.array(list(ao_fracs.values()))[:,1] - wgts[1]/wgts[0]*np.array(list(ao_fracs.values()))[:,0]] #TESTING
    return out

def init_CR_pos(self):
    global CR_top_pos
    CR_top_pos = self.surfaces["sF07"].coefficients["z0"] - self.surfaces["sF04"].coefficients["z0"]
    
def move_CR(self,pos=None):
    global CR_top_pos
    try:
        CR_top_pos
    except NameError:
        init_CR_pos()
    if pos is None:
        pos = CR_top_pos
    # sGRC01 = openmc.ZPlane(6.506+h-half_height) # CR bottom plug 
    # s_CR_pos = openmc.ZPlane(11.365+h-half_height) #For varying CR position
    # Top of self.mats["AIC"] cell is sF07
    l_plug = 11.365-6.506
    pos0 = self.surfaces["sF04"].coefficients["z0"] - l_plug
    posTOP = self.surfaces["sF07"].coefficients["z0"]
    self.surfaces["sGRC01"].coefficients["z0"] = posTOP - pos
    self.surfaces["s_CR_pos"].coefficients["z0"] = posTOP + l_plug - pos
    if openmc.lib.is_initialized:
        openmc.lib.set_surface_z0_by_id(self.surfaces["sGRC01"].id, self.surfaces["sGRC01"].coefficients["z0"])
        openmc.lib.set_surface_z0_by_id(self.surfaces["s_CR_pos"].id, self.surfaces["s_CR_pos"].coefficients["z0"])
    print(f"Moved CR to position {pos - CR_top_pos} cm from fuel top")
    out= {}
    out["materials"] = []
    return out

### For ONE MATERIAL TESTING
# Batch estimated value: 213.6957764742751 +/- 0.48700296567380613
# Moved CR to position 190.51977647427512 cm from fuel top

### Full v3 core
# Batch estimated value: 205.00506784204995 +/- 0.8051954050882285
# Moved CR to position 181.82906784204997 cm from fuel top