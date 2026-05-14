import openmc
import numpy as np
from .core_config import *

def init_geo(self):
    self.cells = {}
    self.universes = {}
    self.fe_pins = []
    self.n_instances = make_core_config(self)
    core_arr = make_core_config2(self)
    # materials.export_to_xml()

    CORE = openmc.RectLattice()
    CORE.lower_left = (-self.fe_pitch*self.n_diam_fe/2,-self.fe_pitch*self.n_diam_fe/2)
    CORE.outer = self.outer
    CORE.pitch = (self.fe_pitch, self.fe_pitch)
    CORE.universes = core_arr
    
    self.cPrism_dist = np.sort(list(dict.fromkeys(self.cPrism_dist)))
    print(self.cPrism_dist)
    cPrism = openmc.model.CruciformPrism(self.cPrism_dist)
    cCORE00a=openmc.Cell(fill=self.mats["cool"], region=-self.surfaces["sF01"])#&+self.surfaces["sF00"])
    cCORE00b=openmc.Cell(fill=self.mats["cool"], region=+self.surfaces["sF08"])#&-self.surfaces["sF09"])
    cCORE01=openmc.Cell(fill=CORE, region=(-cPrism)&+self.surfaces["sF01"]&-self.surfaces["sF08"])
    cCORE02=openmc.Cell(fill=self.mats["HRefl"], region=~(-cPrism)&+self.surfaces["sF01"]&-self.surfaces["sF08"])

    reactor1=openmc.Universe()
    reactor1.add_cells([cCORE00a,cCORE00b,cCORE01,cCORE02])

    cCORE03=openmc.Cell(fill=reactor1, region=-self.surfaces["sCORE02"])
    cCORE04=openmc.Cell(fill=self.mats["ss304"], region=+self.surfaces["sCORE02"]&-self.surfaces["sCORE03"])
    cCORE05=openmc.Cell(fill=self.mats["cool"], region=+self.surfaces["sCORE03"]&-self.surfaces["sCORE04"])
    cCORE06=openmc.Cell(fill=self.mats["SS309L"], region=+self.surfaces["sCORE04"]&-self.surfaces["sCORE05"])
    cCORE07=openmc.Cell(fill=None,region=+self.surfaces["sCORE05"])
    # cCORE_TEST=openmc.Cell(fill=None,region=+sCORE02)

    reactor=openmc.Universe()
    reactor.add_cells([cCORE03,cCORE04,cCORE05,cCORE06,cCORE07])
    # reactor.add_cells([cCORE03,cCORE_TEST]) # HIGH leakage test
    # reactor.add_cells([cCORE03,openmc.Cell(fill=None,region=+sCORE02)])

    cCORE=openmc.Cell(fill=reactor,region=-self.surfaces["sCORE06"])
    #cCORE.plot(pixels=int(1e7), basis="xz",width=(2*(rRPV+100),h+75))
    uni=openmc.Universe()
    uni.add_cells([cCORE])

    geometry = openmc.Geometry()
    geometry.root_universe = uni
    self.geometry = geometry
    # geometry.export_to_xml()