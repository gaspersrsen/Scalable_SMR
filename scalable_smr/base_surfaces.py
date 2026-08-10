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
    self.r_barrel_in = self.r 
    self.r_barrel_out = self.r_barrel_in + self.core_barrel_thickness
    self.h = self.r/0.4699
    self.half_height = 0#(self.h+43.561)/2 #If we want to set center of sim to be in the middle of core

    print(f"Core fuel active height: {self.h}")
    # self.surfaces["sF00"] = openmc.ZPlane(0.0 - self.half_height - 100)
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
    
    ### PLANES FOR FUEL R SECTIONS
    self.r_fuels = []
    for i in range(1,self.n_fuel_r_sections):
        r_n = self.r_fuel * np.sqrt(i / self.n_fuel_r_sections)
        self.r_fuels.append(r_n)
    i_r = 0
    for r in self.r_fuels:
        self.surfaces[f"sFuelR_{i_r}"] = openmc.ZCylinder(0,0,r)
        i_r += 1
    self.r_fuels.append(self.r_fuel)
    self.r_fuels = np.array(self.r_fuels)
    self.surfaces[f"sFuelR_{i_r}"] = openmc.ZCylinder(0,0,self.r_fuel)

    self.surfaces["sF05"] = openmc.ZPlane(self.h_core_top + self.upper_pin_plenum_thickness - self.half_height)
    self.surfaces["sF06"] = openmc.ZPlane(self.h_core_top + self.upper_pin_plenum_thickness + self.upper_pin_cap_thickness - self.half_height)
    self.surfaces["sF07"] = openmc.ZPlane(self.h_core_top + self.upper_pin_plenum_thickness + self.upper_pin_cap_thickness + self.upper_core_cool_thickness - self.half_height)
    self.top_nozzle_pos = self.h_core_top + self.upper_pin_plenum_thickness + self.upper_pin_cap_thickness + self.upper_core_cool_thickness + self.upper_nozzle_thickness
    self.surfaces["sF08"] = openmc.ZPlane(self.top_nozzle_pos - self.half_height)
    self.surfaces["sF09"] = openmc.ZPlane(self.top_nozzle_pos + self.h - self.half_height, boundary_type="vacuum")
    
    self.h_reflector = self.surfaces["sF08"].coefficients["z0"] - self.surfaces["sF01"].coefficients["z0"]

    ### MIXING GRID
    self.mixing_grid_pos = []
    self.n_mix = int((self.surfaces["sF07"].z0-self.surfaces["sF03"].z0)//50) #Between sF02 and sF08
    self.sG=[self.surfaces["sF02"]]
    for i in range(self.n_mix+1):
        self.sG+=[openmc.ZPlane(self.surfaces["sF03"].coefficients["z0"] + i*50 - 5 - self.half_height),
                  openmc.ZPlane(self.surfaces["sF03"].coefficients["z0"] + i*50     - self.half_height)]
        self.mixing_grid_pos += [self.surfaces["sF03"].coefficients["z0"] + i*50 - 5 - self.half_height]
    self.sG+=[self.surfaces["sF08"]]
    # self.rRPV = np.sqrt(1.2/np.pi*(self.n_fe * (self.fe_pitch**2 - np.pi*(25*0.6121**2+264*self.r_zirc**2))
    #                     + (np.pi*self.r**2 - self.n_fe*self.fe_pitch**2)*0.044)
    #                     + (self.r_barrel_in)**2)
    self.rRPV = np.sqrt((self.r_barrel_out)**2 # Outer core barrel
                        + 1.4/np.pi*(self.n_fe*(self.fe_pitch**2 - np.pi*(25*self.surfaces["cyl5"].r**2+264*self.r_zirc**2)) # Fuel element water area
                                    + (np.pi*self.r**2 - self.n_fe*self.fe_pitch**2)*0.044 # Reflector water area
                                     )
                        )
    print(f"Subchannel area: {(self.fe_pitch**2 - np.pi*(25*self.surfaces["cyl5"].r**2+264*self.r_zirc**2))} cm2")
    print(f"RPV inner radius: {self.rRPV}")
    
    self.surfaces["sA01x01"] = openmc.XPlane(10.7518)
    self.surfaces["sA01y01"] = openmc.YPlane(10.7518)
    self.surfaces["sA01x02"] = openmc.XPlane(-10.7518)
    self.surfaces["sA01y02"] = openmc.YPlane(-10.7518)
    
    self.sphere_z0 = self.surfaces["sF01"].z0+np.sin(np.pi/4)*self.rRPV
    self.sphere_inner_r = np.sqrt((np.sin(np.pi/4)*self.rRPV)**2 + self.rRPV**2)
    self.sphere_outer_r = np.sqrt((np.sin(np.pi/4)*self.rRPV)**2 + (self.rRPV+self.RPV_thickness)**2)
    self.sphere_center_z = self.surfaces["sF01"].z0+np.sin(np.pi/4)*self.rRPV
    self.surfaces["ss00"] = openmc.Sphere(z0=self.sphere_center_z, r=self.sphere_inner_r)
    self.surfaces["ss01"] = openmc.Sphere(z0=self.sphere_center_z, r=self.sphere_outer_r, boundary_type="vacuum")
    
    self.surfaces["sCORE_L"] = openmc.YPlane(-self.n_diam_fe*self.fe_pitch)
    self.surfaces["sCORE_R"] = openmc.YPlane(self.n_diam_fe*self.fe_pitch)
    self.surfaces["sCORE_D"] = openmc.XPlane(-self.n_diam_fe*self.fe_pitch)
    self.surfaces["sCORE_U"] = openmc.XPlane(self.n_diam_fe*self.fe_pitch)

    self.surfaces["sCORE02"] = openmc.ZCylinder(0.0, 0.0, self.r_barrel_in)  #INNER CORE BARREL
    self.surfaces["sCORE03"] = openmc.ZCylinder(0.0, 0.0, self.r_barrel_out) #OUTER CORE BARREL
    self.surfaces["sCORE04"] = openmc.ZCylinder(0.0, 0.0, self.rRPV)         #INNER RPV
    self.surfaces["sCORE05"] = openmc.ZCylinder(0.0, 0.0, self.rRPV + self.RPV_thickness, boundary_type="vacuum") #OUTER RPV
    self.surfaces["sCORE06"] = openmc.model.RectangularParallelepiped(-self.rRPV-15,
                                                                      self.rRPV+15,
                                                                      -self.rRPV-15,
                                                                      self.rRPV+15,
                                                                      self.sphere_center_z-self.sphere_outer_r,
                                                                      self.surfaces["sF09"].coefficients["z0"],
                                                    boundary_type="vacuum")


    