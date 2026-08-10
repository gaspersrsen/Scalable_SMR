import openmc
import numpy as np
from .base_surfaces import *
from .mats import *



def make_non_fuel_pins(self):
    make_base_surf(self)
    #Grid_HMP = openmc.model.pin([cyl5,cyl6],[self.mats["cool"],self.mats["inc718"],self.mats["cool"]])
    ###self.surfaces["sF08"]
    Grid_HTP = openmc.model.pin([self.surfaces["cyl5"],self.surfaces["cyl6"]],[self.mats["cool"],self.mats["zirc4"],self.mats["cool"]])
    self.universes["Grid_HTP"] = Grid_HTP

    GRID = openmc.Universe()
    cGRIbot=openmc.Cell(fill=self.mats["NozleBottom"] ,region=-self.surfaces["sF02"])
    cGRItop=openmc.Cell(fill=self.mats["NozleTop"] ,region=+self.surfaces["sF08"])
    
    # ### No mixing grids
    # GRID.add_cells([cGRIbot,cGRItop,
    #                 openmc.Cell(fill=self.mats["cool"] ,region=+self.surfaces["sF02"] &-self.surfaces["sF08"])])
    # ### No mixing grids
    
    cGRIcooltop=openmc.Cell(fill=self.mats["cool"] ,region=+self.sG[-2] &-self.surfaces["sF08"])#
    GRID.add_cells([cGRIbot,cGRItop,cGRIcooltop])
    for i in range(self.n_mix+1):
        cGRIa=openmc.Cell(fill=self.mats["cool"] ,region=+self.sG[2*i]&-self.sG[2*i+1])
        cGRIb=openmc.Cell(fill=self.universes["Grid_HTP"] ,region=+self.sG[2*i+1]&-self.sG[2*i+2])  # Grid 1
        GRID.add_cells([cGRIa, cGRIb])
    
    self.universes["GRID"] = GRID

    self.cells["EndCap"] = openmc.model.pin([self.surfaces["cyl3"]],[self.mats["zirc4"],self.universes["GRID"]])
    self.cells["Fplenum"] = openmc.model.pin([openmc.ZCylinder(0,0,0.0646),self.surfaces["cyl2"],self.surfaces["cyl3"]],[self.mats["inc718"],self.mats["helium"],self.mats["zirc4"],self.universes["GRID"]])
    Gtube = openmc.model.pin([self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])
    self.universes["Gtube"] = Gtube
    Gtube_nozzle = openmc.model.pin([self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["cool"],self.mats["zirc4"],self.mats["NozleTop"]])
    self.universes["Gtube_nozzle"] = Gtube_nozzle
    pinAIC = openmc.model.pin([openmc.ZCylinder(0,0,0.4267),openmc.ZCylinder(0,0,0.4369),openmc.ZCylinder(0,0,0.4839),
                            self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["AIC"],self.mats["helium"],self.mats["ss304"],self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])
    pinAIC_nozzle = openmc.model.pin([openmc.ZCylinder(0,0,0.4267),openmc.ZCylinder(0,0,0.4369),openmc.ZCylinder(0,0,0.4839),
                            self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["AIC"],self.mats["helium"],self.mats["ss304"],self.mats["cool"],self.mats["zirc4"],self.mats["NozleTop"]])
    self.universes["pinAIC"] = pinAIC
    self.universes["pinAIC_nozzle"] = pinAIC_nozzle
    Cplug = openmc.model.pin([openmc.ZCylinder(0,0,0.4839),self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["ss304"],self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])
    self.universes["Cplug"] = Cplug
    # Cplug = openmc.model.pin([cyl4,cyl5],[self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])

    self.cells["cWAT01"]=openmc.Cell(fill=self.mats["NozleBottom"],region=-self.surfaces["sF02"])
    self.cells["cWAT02"]=openmc.Cell(fill=self.mats["cool"],region=+self.surfaces["sF02"]&-self.surfaces["sF07"])
    self.cells["cWAT03"]=openmc.Cell(fill=self.mats["NozleTop"],region=+self.surfaces["sF07"]&-self.surfaces["sF08"])
    self.cells["cWAT04"]=openmc.Cell(fill=self.mats["cool"],region=+self.surfaces["sF08"])

    self.cells["cGTU01"]=openmc.Cell(fill=self.mats["NozleBottom"],region=-self.surfaces["sF02"])
    self.cells["cGTU02"]=openmc.Cell(fill=Gtube,region=+self.surfaces["sF02"]&-self.surfaces["sF07"])
    self.cells["cGTU03"]=openmc.Cell(fill=Gtube_nozzle,region=+self.surfaces["sF07"]&-self.surfaces["sF08"])
    self.cells["cGTU04"]=openmc.Cell(fill=Gtube,region=+self.surfaces["sF08"])#&-self.surfaces["sF09"])


    self.surfaces["sGRC01"] = openmc.ZPlane(self.surfaces["sF04"].z0-4.859) # CR bottom plug 
    self.surfaces["s_CR_pos"] = openmc.ZPlane(self.surfaces["sF04"].z0) #For varying CR position
    self.cells["cGRC01"]=openmc.Cell(fill=self.mats["NozleBottom"],region=-self.surfaces["sF02"])
    self.cells["cGRC02"]=openmc.Cell(fill=Gtube,region=+self.surfaces["sF02"]&-self.surfaces["sGRC01"])
    self.cells["cGRC03"]=openmc.Cell(fill=Cplug,region=+self.surfaces["sGRC01"]&-self.surfaces["s_CR_pos"])
    self.cells["cGRC04"]=openmc.Cell(fill=pinAIC,region=+self.surfaces["s_CR_pos"]&-self.surfaces["sF07"])
    self.cells["cGRC05"]=openmc.Cell(fill=pinAIC_nozzle,region=+self.surfaces["sF07"]&-self.surfaces["sF08"])
    self.cells["cGRC06"]=openmc.Cell(fill=pinAIC,region=+self.surfaces["sF08"])


    WAT = openmc.Universe()
    WAT.add_cells([self.cells["cWAT01"],self.cells["cWAT02"],self.cells["cWAT03"],self.cells["cWAT04"]])
    self.universes["WAT"] = WAT

    GTU = openmc.Universe()
    GTU.add_cells([self.cells["cGTU01"],self.cells["cGTU02"],self.cells["cGTU03"],self.cells["cGTU04"]])
    self.universes["GTU"] = GTU

    GRC = openmc.Universe()
    GRC.add_cells([self.cells["cGRC01"],self.cells["cGRC02"],self.cells["cGRC03"],self.cells["cGRC04"],self.cells["cGRC05"],self.cells["cGRC06"]])
    self.universes["GRC"] = GRC
    

    self.outer = openmc.Universe()
    self.outer.add_cells([openmc.Cell(fill=self.mats["cool"] ,region=-self.surfaces["sCORE04"])])
    self.universes["outer"] = self.outer