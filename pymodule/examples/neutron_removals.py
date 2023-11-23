import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

## define model config object
model_config = {"model":"OLAZR_FM", # model type
                  "projectile":"12c", # projectile, in this case default builtin 12C will be used
                  "target":"12c",
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"classic", # or "none" or "sommerfeld", if not specified none is used
                "level_density":"GC_GEM"
                  }

# now make model from the config
gm = nurex.make_model(model_config)

# we can calculate charge changing cross section for multiple energies
energy = 350.0
cccs = gm.sigma_cc(energy)
incs = gm.sigma_ins(energy)  # incs is array of primary neutron removals, incs[0] corresponds to 1n removal
xncs = gm.sigma_xn(energy)
nevap = gm.n_removals_evaporation() #Total evaporation probability and charged-particle evaporation probabilities for i-neutron removals
s = sum(incs)
print(f"Energy = {energy} MeV/u")
print(f"sigma_cc = {cccs} mb")
print(f"sigma_xn = {xncs} mb")
print(f"sigma_ins = {s} mb")
for i in range(len(incs)):
  Ptot = nevap["Ptot"][i]
  Pch = nevap["Pch"][i]
  Pn = nevap["Pn"][i]
  Pp = nevap["Pp"][i]
  Pd = nevap["Pd"][i]
  Pt = nevap["Pt"][i]
  Phe3 = nevap["Phe3"][i]
  Pa = nevap["Pa"][i]
  Pimf = nevap["Pimf"][i]
  print("------")
  print(f"{i+1}n removal = {incs[i]} mb")
  print(f"Ptot = {Ptot}, Pch = {Pch}, Pn = {Pn}, Pp = {Pp}, Pd = {Pd}, Pt = {Pt}, Phe3 = {Phe3}, Pa = {Pa}, Pimf = {Pimf}")

