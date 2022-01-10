import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

# Example how to scan Fermi Energy coefficient parameter for fensity dependent fermi motion corrections

model_config = {"model":"OLAZR_FM", # Note Constant Fermi Motion coefficients are for NN Cross-section with Fermi motion correction, ie it must be model ending _FM
                "projectile":"28si", 
                "target":"12c",
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"classic", 
                "level_density":"GC_GEM"
                }

energy = 200
print(f"---- Energy = {energy} MeV/u-------")
pm = np.arange(60,100,5) # fermi motion correction range in Mev/c
for c in pm:
  model_config["fermi_motion"] = [c,c]
  gm = nurex.make_model(model_config)
  # we can calculate charge changing cross section for multiple energies  
  rcs = gm.sigma_r(energy)
  cccs = gm.sigma_cc(energy)  
  print(f"FM correction momentum  = {c}, sigma_r = {rcs} mb, sigma_cc = {cccs} mb")  

