import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

# Example how to scan Evaporation Emax parameter

## define model config object
model_config = {"model":"OLAZR_FM",
                "projectile":"48Ca", 
                "target":"12c",
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"classic", 
                "level_density":"GC_GEM"
                }

energy = 400
print(f"---- Energy = {energy} MeV/u-------")

## default evaporation settings
gm = nurex.make_model(model_config)
rcs = gm.sigma_r(energy)
cccs = gm.sigma_cc(energy)  
print(f"Default: sigma_r = {rcs} mb, sigma_cc = {cccs} mb")  

## disable imf evaporation channel
model_config["disable_imf_evaporation"]=1.0
gm = nurex.make_model(model_config)
rcs = gm.sigma_r(energy)
cccs = gm.sigma_cc(energy)  
print(f"Disabled imf: sigma_r = {rcs} mb, sigma_cc = {cccs} mb")  

## disable neuton evaporation channel
model_config["disable_imf_evaporation"]=False
model_config["disable_neutron_evaporation"]=1
gm = nurex.make_model(model_config)
rcs = gm.sigma_r(energy)
cccs = gm.sigma_cc(energy)  
print(f"Dosabled neutron, sigma_r = {rcs} mb, sigma_cc = {cccs} mb")  

