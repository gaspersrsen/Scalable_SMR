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
    upper_pin_plenum_thickness: float = 1.205
    upper_pin_cap_thickness: float = 8.481
    upper_nozzle_thickness: float = 8.481
    upper_core_extra_thickness: float = 9.02
    core_barrel_thickness: float = 5.08
    RPV_thickness: float = 10  # WARNING: 10cm =4in; very low (AP1000 has 8 inches, up to 10'')
    reflector_thickness: float = 18.7174
    initial_boron_ppm: float = 3000.0
    n_diam_fe: int = 7
    one_material: bool = False
    n_fuel_z_sections: int = 6
    without_gadolinia: bool = False
    
    # def __init__(self,
    #     r_fuel=0.4058,
    #     r_helium=0.4140,
    #     r_zirc=0.4750,
    #     pin_pitch=1.2598,
    #     fe_pitch=21.5036,
    #     lower_nozzle_thickness=10.16,
    #     lower_pin_cap_thickness=1.205,
    #     upper_pin_plenum_thickness=1.205,
    #     upper_pin_cap_thickness=8.481,
    #     upper_nozzle_thickness=8.481,
    #     upper_core_extra_thickness=9.02,
    #     core_barrel_thickness=5.08,
    #     RPV_thickness=10,  # WARNING: 10cm =4in; very low (AP1000 has 8 inches, up to 10'')
    #     reflector_thickness=18.7174,
    #     initial_boron_ppm=3000.0,
    #     n_diam_fe=7,
    #     one_material=False,
    #     **kwargs
    # ):
    #     self.parameters:
            
    #     params = locals().copy()
    #     inner_kwargs = params.pop('kwargs')
    #     params.update(inner_kwargs)
        
    def make_model(self):
        import openmc
        import numpy as np
        from .model import init_model
        
        model = init_model(self)
        print(f"Returned a reactor model with: {self.n_diam_fe} fuel elements in diameter.\nEdit, export and run the model.")
        return model