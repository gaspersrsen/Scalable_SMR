import openmc
import numpy as np
from .fuel_elements import *

def get_sim_ij(i,j,d):
    if d%2 == 1:
        if np.abs(i) > np.abs(j):
            return f"{int(np.abs(i)):02d}{int(np.abs(j)):02d}"
        else:
            return f"{int(np.abs(j)):02d}{int(np.abs(i)):02d}"
    else:
        if np.abs(i if i>=0 else i+1) >= np.abs(j if j>=0 else j+1):
            return f"{int(np.abs(i if i>=0 else i+1)):02d}{int(np.abs(j if j>=0 else j+1)):02d}"
        else:
            return f"{int(np.abs(j if j>=0 else j+1)):02d}{int(np.abs(i if i>=0 else i+1)):02d}"



def make_core_config(self):
    init_fuel_elements(self)
    d = self.n_diam_fe
    # global cPrism_dist
    cPrism_dist = [d*self.fe_pitch/2]
    n_instances = []
    list_i = []
    r_core = (d-1) / 2
    core_j = -int(np.floor(d/2))
    for y in np.linspace(-r_core,r_core,d):
        core_i = -int(np.floor(d/2))
        found_edge = False
        for x in np.linspace(-r_core,r_core,d):
            if x >= 0 and y >= 0 and y <= x:
                if filled(x,y,r_core):
                    n_fei = 8
                    if x==0: n_fei/=2
                    if y==0: n_fei/=2
                    if x==y: n_fei/=2
                    n_instances += [int(n_fei)]
                    if self.one_material:
                        continue#pass # ONE MATERIAL
                    elif fatfilled(x,y,r_core) and not self.without_gadolinia:
                        make_new_Gdfuel_element(self, f"{core_i:02d}{core_j:02d}", n_fei,0.0495, 0.08)
                    else: 
                        make_new_fuel_element(self, f"{core_i:02d}{core_j:02d}", int(n_fei), 0.0495)
                else:
                    if not found_edge:
                        if core_i not in list_i:
                            cPrism_dist += [(np.abs(core_j)-.5*float(d%2))*self.fe_pitch]
                            cPrism_dist += [(np.abs(core_i)-.5*float(d%2))*self.fe_pitch]
                            found_edge = True
                            list_i += [core_i]
            core_i += 1
        core_j += 1
    self.cPrism_dist = cPrism_dist
    return n_instances

def make_core_config2(self):
    d = self.n_diam_fe
    if self.one_material:
        make_new_fuel_element(self, f"0000", np.sum(self.n_fe), 0.0495) # ONE MATERIAL
        # self.fes[f"F0000"].plot()
    make_water_element(self)
    core_arr = []
    r_core = (d-1) / 2
    j = -int(np.floor(d/2))
    for y in np.linspace(-r_core,r_core,d):
        row = []
        i = -int(np.floor(d/2))
        for x in np.linspace(-r_core,r_core,d):
            if filled(x,y,r_core):
                if self.one_material:
                    row += [self.fes[f"F0000"]] # ONE MATERIAL
                else:
                    row += [self.fes[f"F{get_sim_ij(i,j,d)}"]]
                
            else:
                row += [self.fes["wWW"]]
            i += 1
        core_arr += [row]
        j += 1
    return core_arr
