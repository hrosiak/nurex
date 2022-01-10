import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

# Example how to scan Fermi Energy coefficient parameter for fensity dependent fermi motion corrections

## define model config object
model_config = {"model":"OLAZR_FMD", # Note Fermi energy coefficient is used only for density dependent Fermi motion correction, ie it must be model ending _FMD
                  "projectile":"28si", 
                  "target":"12c",
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"classic", 
                "level_density":"GC_GEM"
                }

energy = 400
print(f"---- Energy = {energy} MeV/u-------")
fe = np.arange(0.3,0.9,0.02)
for c in fe:
  model_config["fermi_energy_coefficient"] = c
  gm = nurex.make_model(model_config)
  # we can calculate charge changing cross section for multiple energies  
  rcs = gm.sigma_r(energy)
  cccs = gm.sigma_cc(energy)  
  print(f"FE coef = {c}, sigma_r = {rcs} mb, sigma_cc = {cccs} mb")  

