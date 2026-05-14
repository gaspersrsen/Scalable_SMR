import openmc
import numpy as np
from .base_surfaces import *
from .mats import *



def make_non_fuel_pins(self):
    make_base_surf(self)
    #Grid_HMP = openmc.model.pin([cyl5,cyl6],[self.mats["cool"],self.mats["inc718"],self.mats["cool"]])
    ###self.surfaces["sF07"]
    Grid_HTP = openmc.model.pin([self.surfaces["cyl5"],self.surfaces["cyl6"]],[self.mats["cool"],self.mats["zirc4"],self.mats["cool"]])
    self.universes["Grid_HTP"] = Grid_HTP

    GRID = openmc.Universe()
    cGRIbot=openmc.Cell(fill=self.mats["NozleBottom"] ,region=-self.surfaces["sF02"])
    cGRItop=openmc.Cell(fill=self.mats["NozleTop"] ,region=+self.surfaces["sF07"])
    
    ### No mixing grids
    # GRID.add_cells([cGRIbot,cGRItop,
    #                 openmc.Cell(fill=self.mats["cool"] ,region=+self.surfaces["sF02"] &-self.surfaces["sF07"])])
    ### No mixing grids
    
    cGRIcooltop=openmc.Cell(fill=self.mats["cool"] ,region=+self.sG[-2] &-self.surfaces["sF07"])#
    GRID.add_cells([cGRIbot,cGRItop,cGRIcooltop])
    for i in range(self.n_mix+1):
        cGRIa=openmc.Cell(fill=self.mats["cool"] ,region=+self.sG[2*i]&-self.sG[2*i+1])
        cGRIb=openmc.Cell(fill=self.universes["Grid_HTP"] ,region=+self.sG[2*i+1]&-self.sG[2*i+2])  # Grid 1
        GRID.add_cells([cGRIa, cGRIb])
    # GRID.plot(pixels=(2000,2000),color_by='material', colors={self.mats["cool"]:'blue', self.mats["zirc4"]:'cyan', self.mats["inc718"]:'magenta'}, basis='xz')
    # grid_box = openmc.model.RectangularParallelepiped(-self.pin_pitch/2, self.pin_pitch/2, -self.pin_pitch/2, self.pin_pitch/2, self.surfaces["sF02"].z0, self.surfaces["sF07"].z0)
    # cGRIDbox = openmc.Cell(fill=GRID, region=-grid_box)
    # GRID_test = openmc.Universe()
    # GRID_test.add_cells([cGRIDbox])
    self.universes["GRID"] = GRID#_test



    self.cells["EndCap"] = openmc.model.pin([self.surfaces["cyl3"]],[self.mats["zirc4"],self.universes["GRID"]])
    self.cells["Fplenum"] = openmc.model.pin([openmc.ZCylinder(0,0,0.0646),self.surfaces["cyl2"],self.surfaces["cyl3"]],[self.mats["inc718"],self.mats["helium"],self.mats["zirc4"],self.universes["GRID"]])
    Gtube = openmc.model.pin([self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])
    self.universes["Gtube"] = Gtube
    pinAIC = openmc.model.pin([openmc.ZCylinder(0,0,0.4267),openmc.ZCylinder(0,0,0.4369),openmc.ZCylinder(0,0,0.4839),
                            self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["AIC"],self.mats["helium"],self.mats["ss304"],self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])
    self.universes["pinAIC"] = pinAIC
    Cplug = openmc.model.pin([openmc.ZCylinder(0,0,0.4839),self.surfaces["cyl4"],self.surfaces["cyl5"]],[self.mats["ss304"],self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])
    self.universes["Cplug"] = Cplug
    # Cplug = openmc.model.pin([cyl4,cyl5],[self.mats["cool"],self.mats["zirc4"],self.universes["GRID"]])

    self.cells["cWAT01"]=openmc.Cell(fill=self.mats["NozleBottom"],region=-self.surfaces["sF02"])
    self.cells["cWAT02"]=openmc.Cell(fill=self.mats["cool"],region=+self.surfaces["sF02"]&-self.surfaces["sF07"])
    self.cells["cWAT03"]=openmc.Cell(fill=self.mats["NozleTop"],region=+self.surfaces["sF07"])

    self.cells["cGTU01"]=openmc.Cell(fill=self.mats["NozleBottom"],region=-self.surfaces["sF02"])
    self.cells["cGTU02"]=openmc.Cell(fill=Gtube,region=+self.surfaces["sF02"]&-self.surfaces["sF07"])
    self.cells["cGTU03"]=openmc.Cell(fill=self.mats["NozleTop"],region=+self.surfaces["sF07"])


    self.surfaces["sGRC01"] = openmc.ZPlane(6.506+self.h-self.half_height) # CR bottom plug 
    self.surfaces["s_CR_pos"] = openmc.ZPlane(11.365+self.h-self.half_height) #For varying CR position
    self.cells["cGRC01"]=openmc.Cell(fill=self.mats["NozleBottom"],region=-self.surfaces["sF02"])
    self.cells["cGRC02"]=openmc.Cell(fill=Gtube,region=+self.surfaces["sF02"]&-self.surfaces["sGRC01"])
    self.cells["cGRC03"]=openmc.Cell(fill=Cplug,region=+self.surfaces["sGRC01"]&-self.surfaces["s_CR_pos"])
    self.cells["cGRC04"]=openmc.Cell(fill=pinAIC,region=+self.surfaces["s_CR_pos"]&-self.surfaces["sF07"])
    self.cells["cGRC05"]=openmc.Cell(fill=self.mats["NozleTop"],region=+self.surfaces["sF07"])


    WAT = openmc.Universe()
    WAT.add_cells([self.cells["cWAT01"],self.cells["cWAT02"],self.cells["cWAT03"]])
    self.universes["WAT"] = WAT

    GTU = openmc.Universe()
    GTU.add_cells([self.cells["cGTU01"],self.cells["cGTU02"],self.cells["cGTU03"]])
    self.universes["GTU"] = GTU

    GRC = openmc.Universe()
    GRC.add_cells([self.cells["cGRC01"],self.cells["cGRC02"],self.cells["cGRC03"],self.cells["cGRC04"],self.cells["cGRC05"]])
    self.universes["GRC"] = GRC
    

    self.outer = openmc.Universe()
    self.outer.add_cells([openmc.Cell(fill=self.mats["cool"] ,region=-self.surfaces["sCORE04"])])
    self.universes["outer"] = self.outer