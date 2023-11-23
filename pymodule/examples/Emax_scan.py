import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

# Example how to scan Evaporation Emax parameter

## define model config object
model_config = {"model":"OLAZR_FM",
                  "projectile":"28si", 
                  "target":"12c",
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"classic", 
                "level_density":"GC_GEM"
                }

energy = 400
print(f"---- Energy = {energy} MeV/u-------")
fe = np.arange(0,70,5) # scan from 0 in step of 5 MeV
for c in fe:
  model_config["excitation_energy"] = c
  gm = nurex.make_model(model_config)
  rcs = gm.sigma_r(energy)
  cccs = gm.sigma_cc(energy)  
  print(f"FE coef = {c}, sigma_r = {rcs} mb, sigma_cc = {cccs} mb")  

