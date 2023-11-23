import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

# example how to scale neutron removals, ie to match 1n removal cross section

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

# now make model nodel with default setting, ie no scaling
gm = nurex.make_model(model_config)

# we can calculate charge changing cross section for multiple energies
energy = 350.0
cccs = gm.sigma_cc(energy)
incs = gm.sigma_ins(energy)  # incs is array of primary neutron removals, incs[0] corresponds to 1n removal
nevap = gm.n_removals_evaporation() #Total evaporation probability and charged-particle evaporation probabilities for i-neutron removals

print("--- Default Setting ---")
print(f"Energy = {energy} MeV/u")
print(f"sigma_cc = {cccs} mb")
for i in range(len(incs)):
  Pt = nevap["Ptot"][i]  
  nrcs = incs[i] * (1-Pt)
  if(Pt<0):continue
  print(f"{i+1}n removal = {incs[i]} mb, Ptot = {Pt}, {i+1}n removal exp = {nrcs} mb")

n_removal_0 = incs[0] * (1-nevap["Ptot"][0])   # save original neutron removal calculation

## now setup up n-removal scaling
model_config["n_removal_scaling"] = 0.6 # all n-renovals will be scaled by 0.7, choose scaling to reproduce ie 1n removal
gm2 = nurex.make_model(model_config)
cccs = gm2.sigma_cc(energy)
incs = gm2.sigma_ins(energy)  # incs is array of primary neutron removals, incs[0] corresponds to 1n removal
nevap = gm2.n_removals_evaporation() #Total evaporation probability and charged-particle evaporation probabilities for i-neutron removals

print("--- Scaled N-Removals ---")
print(f"Energy = {energy} MeV/u")
print(f"sigma_cc = {cccs} mb")
for i in range(len(incs)):
  Pt = nevap["Ptot"][i]  
  nrcs = incs[i] * (1-Pt)
  if(Pt<0):continue
  print(f"{i+1}n removal = {incs[i]} mb, Ptot = {Pt}, {i+1}n removal exp = {nrcs} mb")


# how to setup scaling to reprodice 1n removal experimental cross section
## not setup up n-removal scaling
n_removal_exp = 42.0  # lets reproduce this
n_scaling = n_removal_exp/n_removal_0
model_config["n_removal_scaling"] = n_scaling # all n-renovals will be scaled by 0.7, choose scaling to reproduce ie 1n removal
gm3 = nurex.make_model(model_config)
cccs = gm3.sigma_cc(energy)
incs = gm3.sigma_ins(energy)  # incs is array of primary neutron removals, incs[0] corresponds to 1n removal
nevap = gm3.n_removals_evaporation() #Total evaporation probability and charged-particle evaporation probabilities for i-neutron removals

print("--- Scaled N-Removals to reproduce 1n removal---")
print(f"n_removal_scaling = {n_scaling}")
print(f"Energy = {energy} MeV/u")
print(f"sigma_cc = {cccs} mb")
for i in range(len(incs)):
  Pt = nevap["Ptot"][i]  
  nrcs = incs[i] * (1-Pt)
  if(Pt<0):continue
  print(f"{i+1}n removal = {incs[i]} mb, Ptot = {Pt}, {i+1}n removal exp = {nrcs} mb")