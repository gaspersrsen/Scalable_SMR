import openmc
import numpy as np
from .core_config import *

def init_geo(self):
    self.cells = {}
    self.fuel_cells_r = []
    self.fuel_cells_r_ids = []
    self.universes = {}
    self.fe_pins = []
    self.fe_instances = {}
    self.fe_enrichments = {}
    self.n_instances = make_core_config(self)
    self.subchannel_list, self.subchannel_connections = generate_octant_subchannels_and_connections(self.n_diam_fe)
    core_arr = make_core_config2(self)

    CORE = openmc.RectLattice()
    CORE.lower_left = (-self.fe_pitch*self.n_diam_fe/2,-self.fe_pitch*self.n_diam_fe/2)
    CORE.outer = self.outer
    CORE.pitch = (self.fe_pitch, self.fe_pitch)
    CORE.universes = core_arr
    self.universes["CORE"] = CORE
    
    self.cPrism_dist = np.sort(list(dict.fromkeys(self.cPrism_dist)))
    print(self.cPrism_dist)
    cPrism = openmc.model.CruciformPrism(self.cPrism_dist)
    self.surfaces["cPrism"] = cPrism
    # cCORE00a = openmc.Cell(fill=self.mats["cool"], region=-self.surfaces["sF01"])
    # self.cells["cCORE00a"] = cCORE00a
    # cCORE00b = openmc.Cell(fill=self.mats["cool"], region=+self.surfaces["sF09"])
    # self.cells["cCORE00b"] = cCORE00b
    cCORE01 = openmc.Cell(fill=CORE, region=(-cPrism)&+self.surfaces["sF01"])
    self.cells["cCORE01"] = cCORE01
    cCORE02 = openmc.Cell(fill=self.mats["HRefl"], region=~(-cPrism)&+self.surfaces["sF01"]&-self.surfaces["sF08"])
    self.cells["cCORE02"] = cCORE02
    cCORE02a = openmc.Cell(fill=self.mats["cool"], region=~(-cPrism)&+self.surfaces["sF08"])
    self.cells["cCORE02a"] = cCORE02a

    reactor1=openmc.Universe()
    # reactor1.add_cells([cCORE00a,cCORE00b,cCORE01,cCORE02])
    reactor1.add_cells([cCORE01,cCORE02,cCORE02a])

    cCORE03 = openmc.Cell(fill=reactor1, region=-self.surfaces["sCORE02"]&+self.surfaces["sF01"])
    self.cells["cCORE03"] = cCORE03
    cCORE04 = openmc.Cell(fill=self.mats["ss304"], region=+self.surfaces["sCORE02"]&-self.surfaces["sCORE03"]&+self.surfaces["sF01"])
    self.cells["cCORE04"] = cCORE04
    cCORE05 = openmc.Cell(fill=self.mats["cool"], region=+self.surfaces["sCORE03"]&-self.surfaces["sCORE04"]&+self.surfaces["sF01"])
    self.cells["cCORE05"] = cCORE05
    cCORE06 = openmc.Cell(fill=self.mats["SS309L"], region=+self.surfaces["sCORE04"]&-self.surfaces["sCORE05"]&+self.surfaces["sF01"])
    self.cells["cCORE06"] = cCORE06
    # cCORE07 = openmc.Cell(fill=None,region=+self.surfaces["sCORE05"])
    # self.cells["cCORE07"] = cCORE07
    # cCORE_TEST=openmc.Cell(fill=None,region=+sCORE02)
    cCORE07 = openmc.Cell(fill=self.mats["cool"], region=-self.surfaces["ss00"]&-self.surfaces["sF01"])
    self.cells["cCORE07"] = cCORE07
    cCORE08 = openmc.Cell(fill=self.mats["SS309L"], region=+self.surfaces["ss00"]&-self.surfaces["ss01"]&-self.surfaces["sF01"])
    self.cells["cCORE08"] = cCORE08

    
    reactor = openmc.Universe()
    self.universes["reactor"] = reactor
    reactor.add_cells([cCORE03,cCORE04,cCORE05,cCORE06,cCORE07,cCORE08])
    # reactor.add_cells([cCORE03,cCORE_TEST]) # HIGH leakage test
    # reactor.add_cells([cCORE03,openmc.Cell(fill=None,region=+sCORE02)])

    cCORE = openmc.Cell(fill=reactor,region=-self.surfaces["sCORE06"])
    self.cells["cCORE"] = cCORE
    #cCORE.plot(pixels=int(1e7), basis="xz",width=(2*(rRPV+100),h+75))
    uni = openmc.Universe()
    self.universes["uni"] = uni
    uni.add_cells([cCORE])

    geometry = openmc.Geometry()
    geometry.root_universe = uni
    self.geometry = geometry
    # geometry.export_to_xml()