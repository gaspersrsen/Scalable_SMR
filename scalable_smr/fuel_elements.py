import openmc
import numpy as np
from .fuel_mat import *
import copy

def init_fuel_elements(self):
    init_fuel_pins(self)
    ### Non fuel cells of pins
    self.cells["cell_1"] = openmc.Cell(fill=self.mats["NozleBottom"], region=-self.surfaces["sF02"])
    self.cells["cell_2"] = openmc.Cell(fill=self.cells["EndCap"], region=+self.surfaces["sF02"]&-self.surfaces["sF03"])
    
    GRID = self.universes["GRID"]
    self.cells["cell_3_He"] = openmc.Cell(fill=self.mats["helium"], region=-self.surfaces["cyl2"]&+self.surfaces["cyl1"]&+self.surfaces["sF03"]&-self.surfaces["sF04"])
    self.cells["cell_3_Zr"] = openmc.Cell(fill=self.mats["zirc4"], region=-self.surfaces["cyl3"]&+self.surfaces["cyl2"]&+self.surfaces["sF03"]&-self.surfaces["sF04"])
    self.cells["cell_3_GRID"] = openmc.Cell(fill=GRID, region=+self.surfaces["cyl3"]&+self.surfaces["sF03"]&-self.surfaces["sF04"])
    
    self.cells["cell_4"] = openmc.Cell(fill=self.cells["Fplenum"], region=+self.surfaces["sF04"]&-self.surfaces["sF05"])
    self.cells["cell_5"] = openmc.Cell(fill=self.cells["EndCap"], region=+self.surfaces["sF05"]&-self.surfaces["sF06"])
    #self.cells["cell_6"] = openmc.Cell(fill=self.mats["cool"], region=+sF06&-sF07)
    self.cells["cell_7"] = openmc.Cell(fill=self.mats["NozleTop"], region=+self.surfaces["sF06"])
    
    self.universes["base_pin"] = openmc.Universe(name="base_pin")
    self.universes["base_pin"].add_cells([self.cells["cell_1"],self.cells["cell_2"],self.cells["cell_3_He"],self.cells["cell_3_Zr"],self.cells["cell_3_GRID"],self.cells["cell_4"],self.cells["cell_5"],self.cells["cell_7"]])


def make_new_fuel_element(self, name, n_instances, enrichment=0.0495):
    print(f"Building fuel element {name} with {n_instances} instances and {enrichment*100} % enrichment")
    fuel_cells = []
    for z_sec in range(self.n_fuel_z_sections):
        new_fuel = make_fuel_material(self, f"{name}_z{z_sec}", enrichment)
        new_fuel.volume = ((17**2-25)*n_instances*self.S_pin*self.z_sec_lengths[z_sec])
        self.mats[f"{name}_z{z_sec}"] = new_fuel
        self.materials += [new_fuel]
        self.fuel_materials += [new_fuel]
        self.Gd_lib[f"{int(new_fuel.id)}"] = 0
        
        fuel_cells += [openmc.Cell(fill=new_fuel, region=-self.surfaces["cyl1"]&+self.surfaces[f"sF03_{z_sec}"]&-self.surfaces[f"sF03_{z_sec+1}"])]

    GRC = self.universes["GRC"]
    GTU = self.universes["GTU"]

    xxx = openmc.Universe(name=f"fuel_pin_{name}")
    xxx.add_cells(fuel_cells)
    xxx.add_cells([openmc.Cell(fill=self.universes["base_pin"], region=+self.surfaces[f"sF00"]&-self.surfaces[f"sF09"])])
    ### NEEDS a new cell, cannot reuse some bug in OpenMC source code
    
    
    new_fel = openmc.RectLattice()
    new_fel.lower_left = (-self.pin_pitch*8.5,-self.pin_pitch*8.5)
    new_fel.outer = self.outer
    new_fel.pitch = (self.pin_pitch,self.pin_pitch)
    new_fel.universes =  [
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, GRC, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, GRC, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GTU, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, GRC, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, GRC, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx]
    ]
    # new_fel.universes = [[xxx for _ in range(17)] for _ in range(17)]
    lc_1=openmc.Cell(fill=new_fel, region=-self.surfaces["sA01x01"]&-self.surfaces["sA01y01"]&+self.surfaces["sA01x02"]&+self.surfaces["sA01y02"])
    # lc_2=openmc.Cell(fill=self.mats["cool"], region=+sA01x01&+sA01y01&-sA01x02&-sA01y02)
    new_fe=openmc.Universe()
    new_fe.add_cells([lc_1])#,lc_2])
    self.fes[f"F{name}"] = new_fe#copy.copy(new_fe)


def make_new_Gdfuel_element(self, name, n_instances, enrichment=0.0495, wGd2O3=0.08):
    Gd_layer_z_stack = []
    fuel_cells = []
    for z_sec in range(self.n_fuel_z_sections):
        new_fuel = make_fuel_material(self, f"{name}_z{z_sec}", enrichment)
        new_fuel.volume = ((17**2-25-16)*n_instances*self.S_pin*self.z_sec_lengths[z_sec])
        self.materials += [new_fuel]
        self.fuel_materials += [new_fuel]
        fuel_cells += [openmc.Cell(fill=new_fuel, region=-self.surfaces["cyl1"]&+self.surfaces[f"sF03_{n_z}"]&-self.surfaces[f"sF03_{n_z+1}"])]
        
        new_fuelGd = make_fuel_material(self, name, enrichment, wGd2O3)
        new_fuelGd.volume = (16*n_instances*self.S_pin*self.z_sec_lengths[z_sec])
        self.fuel_materials += [new_fuelGd]
        self.Gd_lib[f"{int(new_fuelGd.id)}"]=wGd2O3
        Gd_layers = openmc.model.pin([self.surfaces["cyl1"]],[new_fuelGd,None],subdivisions={0:10},divide_vols=True)
        Gd_layer_z_stack += [openmc.Cell(fill=Gd_layers, region=-self.surfaces["cyl1"]&+self.surfaces[f"sF03_{i}"]&-self.surfaces[f"sF03_{i+1}"])]

    GRC = self.universes["GRC"]
    GTU = self.universes["GTU"]
        
    xxx = openmc.Universe(name=f"fuel_pin_{name}")
    xxx.add_cells(fuel_cells)
    xxx.add_cells([openmc.Cell(fill=self.universes["base_pin"], region=+self.surfaces[f"sF00"]&-self.surfaces[f"sF09"])])
    ### NEEDS a new cell, cannot reuse some bug in OpenMC source code    

    FGd = openmc.Universe(name=f"FGd_{name}")
    FGd.add_cells(Gd_layer_z_stack)
    FGd.add_cells([openmc.Cell(fill=self.universes["base_pin"], region=+self.surfaces[f"sF00"]&-self.surfaces[f"sF09"])])
    
    new_fel = openmc.RectLattice() # NEW 16 symetric Gd pins
    new_fel.lower_left = (-self.pin_pitch*8.5,-self.pin_pitch*8.5)
    new_fel.outer = self.outer
    new_fel.pitch = (self.pin_pitch,self.pin_pitch)
    new_fel.universes = [  
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, FGd, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, FGd, xxx, xxx],
         [xxx, xxx, xxx, GRC, xxx, xxx, xxx, FGd, xxx, FGd, xxx, xxx, xxx, GRC, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, FGd, xxx, xxx, xxx, FGd, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, FGd, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, FGd, xxx, xxx, xxx],
         [xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GTU, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx],
         [xxx, xxx, xxx, FGd, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, FGd, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, FGd, xxx, xxx, xxx, FGd, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, GRC, xxx, xxx, xxx, FGd, xxx, FGd, xxx, xxx, xxx, GRC, xxx, xxx, xxx],
         [xxx, xxx, FGd, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, GRC, xxx, xxx, FGd, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],
         [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx]
    ]
    lc_1=openmc.Cell(fill=new_fel, region=-self.surfaces["sA01x01"]&-self.surfaces["sA01y01"]&+self.surfaces["sA01x02"]&+self.surfaces["sA01y02"])
    # lc_2=openmc.Cell(fill=self.mats["cool"], region=+sA01x01&+sA01y01&-sA01x02&-sA01y02)
    new_fe=openmc.Universe()
    new_fe.add_cells([lc_1])#,lc_2])
    self.fes[f"F{name}"] = new_fe#copy.copy(new_fe)


def make_water_element(self):
    # WAT = self.universes["WAT"]
    # wWWL = openmc.RectLattice()
    # wWWL.lower_left = (-self.pin_pitch*9.5,-self.pin_pitch*9.5)
    # wWWL.outer = self.outer
    # wWWL.pitch = (self.pin_pitch,self.pin_pitch)
    # wWWL.universes = [[WAT for i in range(19)]for j in range(19)]
    # L31 = openmc.Cell(fill=wWWL, region=-self.surfaces["sA01x01"]&-self.surfaces["sA01y01"]&+self.surfaces["sA01x02"]&+self.surfaces["sA01y02"])
    # # L32=openmc.Cell(fill=self.mats["cool"], region=+sA01x01&+sA01y01&-sA01x02&-sA01y02)
    # wWW=openmc.Universe()
    # wWW.add_cells([L31])#,L32])
    # self.fes["wWW"] = wWW#copy.copy(wWW)
    
    self.fes["wWW"] = self.outer