from .mats import *
from .geo import *
from .sett import *
from .mixing_functions import *

def init_model(self):
    make_base_mats(self)
    # _ = mat_mix_B(self.initial_boron_ppm) #Initial boron conc
    # if one_material is True:
    #     mat_mix_U(self)
    init_geo(self)
    # make_tallies(self)
    settings_default(self)
    self.model = openmc.model.Model(geometry=self.geometry,materials=self.materials,settings=self.settings)#, tallies=self.tallies)
    return self.model

#TODO move
### CR POSITION
def make_tallies(self):
    self.tallies = openmc.Tallies()
    tallyTest2 = openmc.Tally(tally_id=8889, name="CDI_tally2")
    tallyTest2.scores = ["nu-fission", "absorption", "nu-scatter", "scatter"]
    mats_to_tally = self.cells["cGRC04"].get_all_materials()
    mats_to_tally.update(self.cells["cGRC03"].get_all_materials())
    del mats_to_tally[self.mats["cool"].id]
    tallyTest2.filters = [openmc.CellFilter([self.cells["cGRC04"],self.cells["cGRC03"]]), openmc.MaterialFilter(list(mats_to_tally))]
    self.tallies += [tallyTest2]