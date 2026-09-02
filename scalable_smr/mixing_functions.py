from .geo import *
import openmc.data
import numpy as np

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


# ### BORON MIXING
# def whole_boron_mix(self, conc):
#     global mix_ao
    
#     cool_positions = []
#     for _pos, mat in enumerate(self.model.materials):
#         if mat.name == "cool":
#             cool_positions.append(_pos)
#     cool = self.mats["cool"]
#     pure_water = self.mats["pure_water"]
#     pure_boric_acid = self.mats["pure_boric_acid"]
#     pure_water_density = pure_water.get_mass_density()
#     pure_boric_acid_density = pure_boric_acid.get_mass_density()
    
#     # Convert ppm back to raw mass fraction for the scaling formulas
#     C_mass = conc * 1e-6
    
#     mix_ao, ao_fracs, wgts_mix = openmc.search.get_ao_mix_materials(
#         [pure_water, pure_boric_acid], [None, C_mass],
#         fracs_target=[None, "B"],
#         percent_type='wo', return_wgts=True)
    

#     new_mats = []
#     for mat in [cool]:
#         mat.update_material(mix_ao,
#                             rho=pure_water_density*wgts_mix[0]+pure_boric_acid_density*wgts_mix[1],
#                             rho_units='g/cm3'
#                             )
#         new_mats.append(mat)
#     rho_mix = cool.get_mass_density()
#     m_boron = np.sum([cool.get_mass_density(nuc) for nuc in ["B10","B11"]])
#     print(f"Rebuilt material with {m_boron/rho_mix*1e6} ppm")
    
#     for _pos, cool_pos in enumerate(cool_positions):
#         self.model.materials[cool_pos] = new_mats[_pos]
    
#     out = {}
#     out["materials"] = [cool]
    
#     # 1. Safely extract the exact nuclide ordering used by the active CDI tally
#     tally_nucs = []
#     for t in self.model.tallies:
#         if t.id == 8889:
#             tally_nucs = t.nuclides
#             break
            
#     if not tally_nucs:
#         tally_nucs = list(mix_ao.keys())
    
#     nuc_dens_water = pure_water.get_nuclide_atom_densities()
#     nuc_dens_boric = pure_boric_acid.get_nuclide_atom_densities()
    
#     # Extract the true physical scalar volume fraction of boric acid (index 1)
#     v_b = wgts_mix[1]
    
#     # 2. Compute the exact Non-Linear Volumetric Correction Factor (Chain Rule link)
#     density_ratio = pure_water_density / pure_boric_acid_density
#     mass_to_vol_correction = 1.0 / (1.0 + C_mass * (density_ratio - 1.0))
    
#     wgts_out = []
#     # 3. Key-matching loop applies the correction to ensure perfect tally alignment
#     for nuc in tally_nucs:
#         N_mix = mix_ao.get(nuc, 0.0)
        
#         # Pull unmixed constituent baseline atom densities
#         N_b = nuc_dens_boric.get(nuc, 0.0)
#         N_w = nuc_dens_water.get(nuc, 0.0)
        
#         # Volumetric sensitivity base (v_b is now a clean scalar)
#         vol_sensitivity = (v_b * (N_b - N_w)) / N_mix if N_mix > 0 else 0.0
        
#         # Package the non-linear transformation directly into the flat output weight
#         final_mass_sensitivity = vol_sensitivity * mass_to_vol_correction
#         wgts_out.append(final_mass_sensitivity)
        
#     wgts_out = np.array(wgts_out)
#     # print("Packaged Non-Linear Mass Sensitivity Weights:\n", wgts_out)
    
#     out["nuc_fractions"] = [wgts_out]
#     return out


def ideal_boron_mix(self, conc):
    """
    Idealized phase-blending function computing nuc_fractions as the true
    logarithmic sensitivity derivative (accounting for water displacement).
    """
    M_B = openmc.data.atomic_weight("B")
    M_H3BO3 = self.mats["pure_boric_acid"].average_molar_mass * 7
    rho_boric_solid = 1.435
    
    cool = self.mats["cool"]
    pure_water = self.mats["pure_water"]
    pure_water_density = pure_water.get_mass_density()
    
    C_mass_fraction = conc * 1e-6
    w_boric_mass = C_mass_fraction * (M_H3BO3 / M_B)
    w_water_mass = 1.0 - w_boric_mass
    
    rho_mix_ideal = 1.0 / ((w_water_mass / pure_water_density) + (w_boric_mass / rho_boric_solid))
    
    # Core Material Blending
    mix_ao, ao_fracs, wgts_mix = openmc.search.get_ao_mix_materials(
        [pure_water, self.mats["pure_boric_acid"]], 
        [None, C_mass_fraction],
        fracs_target=[None, "B"],
        percent_type='wo', return_wgts=True
    )
    
    cool.update_material(mix_ao, rho=rho_mix_ideal, rho_units='g/cm3')
    
    cool_positions = [i for i, m in enumerate(self.model.materials) if m.name == "cool"]
    for pos in cool_positions:
        self.model.materials[pos] = cool

    # Tally Alignment Track
    tally_nucs = []
    for t in self.model.tallies:
        if t.id == 8889:
            tally_nucs = t.nuclides
            break
    if not tally_nucs:
        tally_nucs = list(mix_ao.keys())
        
    # Extract structural phase densities
    nuc_dens_water = pure_water.get_nuclide_atom_densities()
    nuc_dens_boric = self.mats["pure_boric_acid"].get_nuclide_atom_densities()
    
    # wgts_mix[1] represents the precise volume fraction of boric acid (v_b)
    v_b = wgts_mix[1]
    
    # Analytical mass-to-volume non-linear scaling correction derivative
    # d(ln v_b) / d(ln C) = 1 - v_b * (1 - rho_boric / rho_water)
    mass_to_vol_correction = 1.0 - v_b * (1.0 - (rho_boric_solid / pure_water_density))
    
    wgts_out = []
    for nuc in tally_nucs:
        N_mix = mix_ao.get(nuc, 0.0)
        N_b = nuc_dens_boric.get(nuc, 0.0)
        N_w = nuc_dens_water.get(nuc, 0.0)
        
        # Volumetric derivative component: v_b * (N_boric - N_water) / N_mix
        vol_sensitivity = (v_b * (N_b - N_w)) / N_mix if N_mix > 0 else 0.0
        
        # Total chain-rule coupled sensitivity derivative
        final_sensitivity = vol_sensitivity * mass_to_vol_correction
        wgts_out.append(final_sensitivity)
        
    print(f"Ideal Rebuild: {conc:.1f} ppm | Ideal Density: {rho_mix_ideal:.5f} g/cc")
    return {"materials": [cool], "nuc_fractions": [np.array(wgts_out, dtype=float)]}


def paper_backed_boron_mix(self, conc, temperature_k=573.15):
    """
    Advanced non-ideal boron mixing function computing nuc_fractions as a flat 1D array
    representing the exact atomic source contribution from boric acid (w_boric).
    """
    M_B = openmc.data.atomic_weight("B")
    M_H3BO3 = self.mats["pure_boric_acid"].average_molar_mass * 7
    
    cool = self.mats["cool"]
    pure_water = self.mats["pure_water"]
    pure_water_density = pure_water.get_mass_density()
    
    # Temperature-Dependent Apparent Molar Volume Fit (Abdulagatov 2004)
    T = temperature_k
    V_phi = (-208.675 
             + 1.83984 * T 
             - 4.24925e-3 * T**2 
             + 3.32831e-6 * T**3)
    
    def compute_aqueous_density(target_ppm):
        w_B = target_ppm * 1e-6
        moles_B = w_B / M_B
        mass_water = 1.0 - (w_B * (M_H3BO3 / M_B))
        m = moles_B / mass_water if mass_water > 0 else 0.0
        
        numerator = 1000.0 + m * M_H3BO3
        denominator = (1000.0 / pure_water_density) + m * V_phi
        return numerator / denominator / 1000.0
    
    C_mass_fraction = conc * 1e-6
    mix_ao, ao_fracs, wgts_mix = openmc.search.get_ao_mix_materials(
        [pure_water, self.mats["pure_boric_acid"]], 
        [None, C_mass_fraction],
        fracs_target=[None, "B"],
        percent_type='wo', return_wgts=True
    )
    
    rho_mix = compute_aqueous_density(conc)
    cool.update_material(mix_ao, rho=rho_mix, rho_units='g/cm3')
    
    cool_positions = [i for i, m in enumerate(self.model.materials) if m.name == "cool"]
    for pos in cool_positions:
        self.model.materials[pos] = cool

    # Tally Alignment Track
    tally_nucs = []
    for t in self.model.tallies:
        if t.id == 8889:
            tally_nucs = t.nuclides
            break
    if not tally_nucs:
        tally_nucs = list(mix_ao.keys())
        
    # Generate 1D Boric Acid Fraction Split [w_boric]
    nuc_dens_water = pure_water.get_nuclide_atom_densities()
    nuc_dens_boric = self.mats["pure_boric_acid"].get_nuclide_atom_densities()
    
    wgts_out = []
    for nuc in tally_nucs:
        atoms_from_water = wgts_mix[0] * nuc_dens_water.get(nuc, 0.0)
        atoms_from_boric = wgts_mix[1] * nuc_dens_boric.get(nuc, 0.0)
        total_atoms = atoms_from_water + atoms_from_boric
        
        w_boric = atoms_from_boric / total_atoms if total_atoms > 0 else 0.0
        wgts_out.append(w_boric)
        
    print(f"Non-Ideal Rebuild: {conc:.1f} ppm | Sol Density: {rho_mix:.5f} g/cc | V_phi: {V_phi:.2f} cm3/mol")
    return {"materials": [cool], "nuc_fractions": [np.array(wgts_out, dtype=float)]}


def init_CR_pos(self):
    global CR_top_pos
    CR_top_pos = self.surfaces["sF08"].coefficients["z0"] - self.surfaces["sF04"].coefficients["z0"]
    
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
    l_plug = 11.365-6.506
    pos0 = self.surfaces["sF04"].coefficients["z0"] - l_plug
    posTOP = self.surfaces["sF08"].coefficients["z0"]
    self.surfaces["sGRC01"].coefficients["z0"] = posTOP - pos
    self.surfaces["s_CR_pos"].coefficients["z0"] = posTOP + l_plug - pos
    if openmc.lib.is_initialized:
        openmc.lib.set_surface_z0_by_id(self.surfaces["sGRC01"].id, self.surfaces["sGRC01"].coefficients["z0"])
        openmc.lib.set_surface_z0_by_id(self.surfaces["s_CR_pos"].id, self.surfaces["s_CR_pos"].coefficients["z0"])
    print(f"Moved CR to position {pos - CR_top_pos} cm from fuel top")
    out= {}
    out["materials"] = []
    return out