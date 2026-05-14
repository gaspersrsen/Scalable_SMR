import openmc

def settings_default(self):
   settings = openmc.Settings()
   settings.batches = 100
   # settings.inactive = 20
   settings.particles = 10000#n_particles
   # settings.generations_per_batch = 1#generations_per_batch
   settings.ptables = True
   # settings.verbosity = 4
   settings.temperature = {"method":"interpolation","tolerance":100}
   settings.output = {"tallies":False}
   # settings.source = openmc.IndependentSource(constraints = {'fissionable': True})
   # settings.source = openmc.IndependentSource(space = openmc.stats.Box(lower_left = geometry.bounding_box.lower_left,
   #                                                                   upper_right = geometry.bounding_box.upper_right,
   #                                                                   only_fissionable = True),
   #                                            constraints = {'fissionable': True})
   settings.source = openmc.IndependentSource(
                                          space = openmc.stats.Box(lower_left = (-self.fe_pitch*self.n_diam_fe/2, -self.fe_pitch*self.n_diam_fe/2, -self.half_height),
                                                                  upper_right = (self.fe_pitch*self.n_diam_fe/2,  self.fe_pitch*self.n_diam_fe/2,  self.half_height)),
                                          constraints = {'fissionable': True})
   settings.export_to_xml()
   self.settings = settings