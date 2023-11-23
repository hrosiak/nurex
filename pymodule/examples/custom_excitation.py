import matplotlib.pyplot as plt
import numpy as np
import pynurex as nurex

## This is example how to use custom excitation function
## Custom excitation function is specified for each i-hole state
## It is supplied using array of excitation energies and array of corresponding distribution value
## Provided values are normalized to 1.0 and fitted using cubic splines 
## when not specified custom excitation function is 0.0

model_config = {"model":"OLAZR", # model type
                        "projectile":"12c", # projectile, in this case default builtin 12C will be used
                        "target":"12c",
                        "charge_changing_correction":"none",
                        "coulomb_correction":"classic", # or "none" or "sommerfeld", if not specified none is used
                        "level_density":"GC_GEM"
                        }

energy = 1000.
# get sigma cc without evaporation                        
gm = nurex.make_model(model_config)
sigma_noevap = gm.sigma_cc(energy)

# get sigma cc with default evaporation, ie GS model
model_config["charge_changing_correction"] = "evaporation"
gm = nurex.make_model(model_config)
sigma_gs = gm.sigma_cc(energy)


## now switch to custom excitation function 
## by default without setting it up it returns 0.0
## ie no evaporation probability
gm.set_excitation_function_type(nurex.excitation_function_type.CUSTOM)  # to switch back to default use nurex.excitation_function_type.GS
sigma_custom0 = gm.sigma_cc(1000)

print(f"Energy = {energy} MeV/u")
print(f"CC sigma_noevap = {sigma_noevap}, sigma_gs={sigma_gs}, sigma_custom0 = {sigma_custom0}")

## now setting up custom excitation
## we need to supply array of energies and corresponding probabilities, ie from data file or other model
## next we just fill from GS model
energies = np.linspace(0,80,400) # array of 400 energies from 0 to 80 MeV
nurex.custom_excitation_reset()
for i in range(1,5):   # loop over holes from 1 to 4
    ex = nurex.Emax(gm.projectile(),gm.get_evaporation_parameters()) # this is default Emax parameter used in nurex
    w = [nurex.w_gs(e,i,ex) for e in energies]  # generate probabilities for energies array
    nurex.custom_excitation_set(energies=energies, w=w, i=i)  # now setting up for i-th hole

sigma_customGS = gm.sigma_cc(1000)
print(f"CC sigma_customGS = {sigma_customGS}")

## now squeezing GS distribution into half exc. energy range
## this shift exc. energies distribution to lower energies 
## therefore supress evaporation probabilities, ie eveporation correction will be smaller
nurex.custom_excitation_reset() # this is importtant to reset previous setup in case we are modifying it
energies = np.linspace(0,40,400)
for i in range(1,5):   
    w = [nurex.w_gs(2*e,i,ex) for e in energies]
    nurex.custom_excitation_set(energies=energies, w=w, i=i)
sigma_custom = gm.sigma_cc(1000)
print(f"CC sigma_custom = {sigma_custom}")

## this is example how to print and cross check 
energies = np.linspace(0,50,500)
w1 = [nurex.custom_excitation_w(energy=e, i=1) for e in energies]
w2 = [nurex.custom_excitation_w(energy=e, i=2) for e in energies]
w3 = [nurex.custom_excitation_w(energy=e, i=3) for e in energies]
w4 = [nurex.custom_excitation_w(energy=e, i=4) for e in energies]
plt.figure()
plt.title("Excitation energy function")
plt.plot(energies, w1, label="1 hole")
plt.plot(energies, w2, label="2 hole")
plt.plot(energies, w3, label="3 hole")
plt.plot(energies, w4, label="4 hole")
plt.xlabel("Excitation energy [MeV]")
plt.grid()
plt.legend()
plt.show()