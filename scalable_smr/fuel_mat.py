import openmc
import numpy as np
from .non_fuel_pins import *

def init_fuel_pins(self):
    make_non_fuel_pins(self)
    self.fes = {} #Fuel element dictionary
    self.fuel_materials = openmc.Materials()
    self.S_pin = (0.4058**2*np.pi)
    self.Gd_lib={}

def get_ao_fuel(w, wGd=0.0):
    ### Look at mole fraction: https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)#Mole_fraction
    mO16 = openmc.data.atomic_mass("O16")
    mU235 = openmc.data.atomic_mass("U235")
    mU238 = openmc.data.atomic_mass("U238")
    mGd_array = np.array([151.9198,153.9209,154.9226,155.9221,156.923,157.9241,159.9271])
    Gd_abundance = np.array([0.0020,0.0218,0.1480,0.2047,0.1565,0.2484,0.2186])
    mGd = np.sum(mGd_array*Gd_abundance) 
    AVOGADRO = 6.0221E23
    BARN = 1e-24
    rho = (10.962 - 3.48*wGd) * 0.97 #Theoretical dens * sintering eff, valid for <12%w/o Gd2O3
    
    if wGd==0.0:
        xU235 = w/mU235 / (w/mU235 + (1-w)/mU238)
        xU238 = 1 - xU235
        mU = xU235*mU235 + xU238*mU238
        mUO2 = mU + 2.0 * mO16
        ao_tot = (rho / mUO2) * AVOGADRO #number of MOLECULES / volume
        UO2_fuel = np.zeros(3)
        UO2_fuel[0] = ao_tot * xU235 * BARN   #U-235
        UO2_fuel[1] = ao_tot * xU238 * BARN   #U-238
        UO2_fuel[2] = 2.0 * ao_tot * BARN     # O-16
        N_tot = np.sum(UO2_fuel)
        return [UO2_fuel, N_tot]
    else:
        #Atomic fractions in UO2
        xU235 = w/mU235 / (w/mU235 + (1-w)/mU238)
        xU238 = 1 - xU235
        #Molar masses
        mU = xU235*mU235 + xU238*mU238
        mGd2O3 = 2.0 * mGd + 3.0 * mO16
        mUO2 = mU + 2.0 * mO16
        #Atomic fractions in Mixture
        wGd2O3 = wGd/mGd / (wGd/mGd + (1-wGd)/mUO2)
        xGd2O3 = wGd2O3/mGd2O3 / (wGd2O3/mGd2O3 + (1-wGd2O3)/mUO2)
        xU = 1 - xGd2O3
        #Molar mass of mixture
        mMix = mUO2 * xU + mGd2O3 * xGd2O3
        #Total atom density
        ao_tot = rho * AVOGADRO / mMix
        #atom number concentration
        Gd_fuel = np.zeros(10)
        Gd_fuel[0] = ao_tot * xU * xU235 * BARN  #U-235
        Gd_fuel[1] = ao_tot * xU * xU238 * BARN  #U-238
        ao_O = (2.0 * ao_tot * xU) + (3.0/2 * ao_tot * xGd2O3) # O-16
        Gd_fuel[2] = ao_O * BARN    # O-16
        Gd_fuel[3:] = np.array(ao_tot / 2 * xGd2O3 * Gd_abundance * BARN)
        N_tot = np.sum(Gd_fuel)
        return [Gd_fuel, N_tot]
    
def make_fuel_material(self, name, w, wGd2O3=0.0):
    if wGd2O3==0.0:
        [UO2_fuel, N_tot] = get_ao_fuel(w, wGd2O3)
        new_fuel = openmc.Material(name=f"fuel{name}", temperature=900)
        new_fuel.add_nuclide("U235",UO2_fuel[0],"ao")
        new_fuel.add_nuclide("U238",UO2_fuel[1],"ao")
        new_fuel.add_nuclide("O16",UO2_fuel[2],"ao")
        new_fuel.depletable=True
        new_fuel.set_density("atom/b-cm", N_tot)
        return new_fuel
    else:
        [Gd_fuel, N_tot] = get_ao_fuel(w, wGd2O3)
        new_fuelGd = openmc.Material(name=f"fuelGd{name}", temperature=900)
        new_fuelGd.add_nuclide("U235",Gd_fuel[0],"ao")
        new_fuelGd.add_nuclide("U238",Gd_fuel[1],"ao")
        new_fuelGd.add_nuclide("O16",Gd_fuel[2],"ao")
        new_fuelGd.add_nuclide("Gd152", Gd_fuel[3],"ao")
        new_fuelGd.add_nuclide("Gd154", Gd_fuel[4],"ao")
        new_fuelGd.add_nuclide("Gd155", Gd_fuel[5],"ao")
        new_fuelGd.add_nuclide("Gd156", Gd_fuel[6],"ao")
        new_fuelGd.add_nuclide("Gd157", Gd_fuel[7],"ao")
        new_fuelGd.add_nuclide("Gd158", Gd_fuel[8],"ao")
        new_fuelGd.add_nuclide("Gd160", Gd_fuel[9],"ao")
        new_fuelGd.depletable=True
        new_fuelGd.set_density("atom/b-cm", N_tot)
        return new_fuelGd