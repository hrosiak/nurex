import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

## define model config object
model_config = {"model":"OLAZR_FM", # model type
                  "projectile":"12c", # projectile, in this case default builtin 12C will be used
                  "target":{          # target is also 12C but now with custom nuclear densitiesm in this case HO type
                  "nucleus":"12c",
                  "proton_density":{"type":"ho", "parameters":[1.548415436,1.6038565]},
                  "neutron_density":{"type":"ho", "parameters":[1.548415436,1.6038565]}
                  },
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"classic", # or "none" or "sommerfeld", if not specified none is used
                "level_density":"GC_GEM"
                  }

# now make model from the config
for proj in ["13c","12c","12c","13c"]:
  print("-----------")
  model_config["projectile"] = proj
  gm = nurex.make_model(model_config)
  # we can calculate charge changing cross section for multiple energies
  energy = 150
  rcs = gm.sigma_r(energy)
  cccs = gm.sigma_cc(energy)
  incs = gm.sigma_ins(energy)  # incs is array of primary neutron removals, incs[0] corresponds to 1n removal
  nevap = gm.n_removals_evaporation() #Total evaporation probability and charged-particle evaporation probabilities for i-neutron removals

  print(f"Energy = {energy} MeV/u")
  print(f"sigma_r = {rcs} mb")
  print(f"sigma_cc = {cccs} mb")
  for i in range(len(incs)):
    Pt = nevap["Ptot"][i]
    Pch = nevap["Pch"][i]
    print(f"{i+1}n removal = {incs[i]} mb, Ptot = {Pt}, Pch = {Pch}")

