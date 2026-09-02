from dataclasses import dataclass

@dataclass
class SMR:
    r_fuel: float = 0.4058
    r_helium: float = 0.4140
    r_zirc: float = 0.4750
    pin_pitch: float = 1.2598
    fe_pitch: float = 21.5036
    lower_nozzle_thickness: float = 10.16
    lower_pin_cap_thickness: float = 1.205
    upper_pin_plenum_thickness: float = 13.490
    upper_pin_cap_thickness: float = 1.205
    upper_nozzle_thickness: float = 9.020
    upper_core_cool_thickness: float = 8.481
    reflector_thickness: float = 18.7174
    core_barrel_thickness: float = 5.08 # 2 in
    RPV_thickness: float = 10.16  # WARNING: 10cm =4in; very low (AP1000 has 8 inches, up to 10'')
    initial_boron_ppm: float = 3000.0
    n_diam_fe: int = 7
    one_material: bool = False
    n_fuel_z_sections: int = 6
    n_fuel_r_sections: int = 1
    without_gadolinia: bool = False
    fuel_enrichment: float = 0.0495
    enrichment_function: callable = None
    single_r_fuel: bool = False
        
    def make_model(self):
        import openmc
        import numpy as np
        from .model import init_model
        
        model = init_model(self)
        print(f"Returned a reactor model with: {self.n_diam_fe} fuel elements in diameter.\nEdit, export and run the model.")
        return model
    
    def set_boron(self, conc):
        from .mixing_functions import ideal_boron_mix
        out = ideal_boron_mix(self, conc)
        return out