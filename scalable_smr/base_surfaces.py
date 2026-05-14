import openmc
import numpy as np
from .fe_instance_counter import *


def make_base_surf(self):
    self.n_fe = get_n_fe(self.n_diam_fe) #Number of fuel elements in the core, used for calculating the RPV radius
    self.surfaces = {}

    self.surfaces["cyl1"] = openmc.ZCylinder(0,0,self.r_fuel)
    self.surfaces["cyl2"] = openmc.ZCylinder(0,0,self.r_helium)
    self.surfaces["cyl3"] = openmc.ZCylinder(0,0,self.r_zirc)
    self.surfaces["cyl4"] = openmc.ZCylinder(0,0,0.5715)
    self.surfaces["cyl5"] = openmc.ZCylinder(0,0,0.6121)
    self.surfaces["cyl6"] = openmc.ZCylinder(0,0,0.6347)
    

    self.r = self.reflector_thickness + self.n_diam_fe*self.fe_pitch/2
    self.dr = 99.060 - 93.98
    self.h = self.r/0.4699
    self.half_height = (self.h+43.561)/2

    print(f"Core fuel active height: {self.h}")
    self.surfaces["sF00"] = openmc.ZPlane(0.0 - self.half_height - self.h/2)
    self.surfaces["sF01"] = openmc.ZPlane(0.0 - self.half_height)
    self.surfaces["sF02"] = openmc.ZPlane(self.lower_nozzle_thickness - self.half_height)
    self.h_core_bot = self.lower_nozzle_thickness + self.lower_pin_cap_thickness
    self.surfaces["sF03"] = openmc.ZPlane(self.h_core_bot - self.half_height)
    self.h_core_top = self.h_core_bot + self.h
    self.surfaces["sF04"] = openmc.ZPlane(self.h_core_top - self.half_height)
    
    ### PLANES FOR FUEL Z SECTIONS
    def flux_height(h):
        if h < 0:
            return 0
        elif h < self.h:
            return 0.3+0.7*np.sin(np.pi*(h)/self.h)
        else:
            return 0
        
    sec_wgts = []
    for x in range(1,self.n_fuel_z_sections+1):
        sec_wgts.append(flux_height(x*self.h/(self.n_fuel_z_sections+1))**2)
        
    sec_wgts = np.array(sec_wgts)
    self.z_sec_lengths = sec_wgts/np.sum(sec_wgts)*self.h
    z_sec_surf_pos = np.cumsum(self.z_sec_lengths)
    print(f"Z section surface positions: {z_sec_surf_pos}")
    for i in range(self.n_fuel_z_sections):
        self.surfaces[f"sF03_{i+1}"] = openmc.ZPlane(self.h_core_bot + z_sec_surf_pos[i]- self.half_height)
        # self.surfaces[f"sF03_{i}"] = openmc.ZPlane(z_sec_surf_pos[i-1]) #These are used in local coordinates
    self.surfaces[f"sF03_0"] = self.surfaces["sF03"]
    self.surfaces[f"sF03_{self.n_fuel_z_sections+1}"] = self.surfaces["sF04"]
    ### PLANES FOR FUEL Z SECTIONS
    
    self.surfaces["sF05"] = openmc.ZPlane(self.h_core_top + self.upper_pin_plenum_thickness - self.half_height)
    self.surfaces["sF06"] = openmc.ZPlane(self.h_core_top + self.upper_pin_plenum_thickness + self.upper_pin_cap_thickness - self.half_height)
    self.top_nozzle_pos = self.h_core_top + self.upper_pin_plenum_thickness + self.upper_pin_cap_thickness + self.upper_nozzle_thickness
    self.surfaces["sF07"] = openmc.ZPlane(self.top_nozzle_pos - self.half_height)
    self.surfaces["sF08"] = openmc.ZPlane(self.top_nozzle_pos + self.upper_core_extra_thickness - self.half_height)
    self.surfaces["sF09"] = openmc.ZPlane(self.top_nozzle_pos + self.upper_core_extra_thickness - self.half_height + self.h/2)

    ### MIXING GRID
    self.n_mix = int((self.surfaces["sF08"].z0-self.surfaces["sF02"].z0)//50) #Between sF02 and sF07
    self.sG=[self.surfaces["sF02"]]
    for i in range(self.n_mix+1):
        self.sG+=[openmc.ZPlane(self.lower_nozzle_thickness + i*50 - 5 - self.half_height),
                  openmc.ZPlane(self.lower_nozzle_thickness + i*50     - self.half_height)]
    self.sG+=[self.surfaces["sF07"]]
    self.rRPV = np.sqrt(1.2/np.pi*(self.n_fe * (self.fe_pitch**2 - np.pi*(25*0.6121**2+264*self.r_zirc**2))
                        + (np.pi*self.r**2 - self.n_fe*self.fe_pitch**2)*0.044)
                        + (self.r+self.dr)**2)
    print(f"RPV inner radius: {self.rRPV}")
    self.surfaces["sCORE_L"] = openmc.YPlane(-self.n_diam_fe*self.fe_pitch)
    self.surfaces["sCORE_R"] = openmc.YPlane(self.n_diam_fe*self.fe_pitch)
    self.surfaces["sCORE_D"] = openmc.XPlane(-self.n_diam_fe*self.fe_pitch)
    self.surfaces["sCORE_U"] = openmc.XPlane(self.n_diam_fe*self.fe_pitch)

    self.surfaces["sCORE02"] = openmc.ZCylinder(0.0, 0.0, self.r) 
    self.surfaces["sCORE03"] = openmc.ZCylinder(0.0, 0.0, self.r + self.reflector_thickness) #OUTER CORE BARREL
    self.surfaces["sCORE04"] = openmc.ZCylinder(0.0, 0.0, self.rRPV)
    self.surfaces["sCORE05"] = openmc.ZCylinder(0.0, 0.0, self.rRPV + self.RPV_thickness) #OUTER RPV
    self.surfaces["sCORE06"] = openmc.model.RectangularParallelepiped(-self.rRPV-15, self.rRPV+15, -self.rRPV-15, self.rRPV+15, 0.0-self.half_height-self.h/2, 43.561+self.h-self.half_height+self.h/2,
                                                    boundary_type="vacuum")

    self.surfaces["sA01x01"] = openmc.XPlane(10.7518)
    self.surfaces["sA01y01"] = openmc.YPlane(10.7518)
    self.surfaces["sA01x02"] = openmc.XPlane(-10.7518)
    self.surfaces["sA01y02"] = openmc.YPlane(-10.7518)
    
    