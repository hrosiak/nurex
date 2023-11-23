import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np
import pprint
pp = pprint.PrettyPrinter(indent=4).pprint

# get predefined 12C nucleus densities
# nucleus contains separate density for protons and neutrons
nuc = nurex.default_nucleus("12c")
dproton = nuc.get_density_proton()
dneutron = nuc.get_density_proton()

# check Rrms of protons and neutrons
Rp = nurex.r_rms(dproton)
Rn = nurex.r_rms(dproton)
print("Rrms proton = ", Rp)
print("Rrms neutron = ", Rn)

# get config dictionary for proton ednsity function
# this can be used to make new density
config_proton = nurex.density_object(dproton)
print("Config object for proton density:")
pp(config_proton)

# get config object for whole nucleus, this can be used for Glauber model config object
config_nucleus = nurex.nucleus_object(nuc)
print("Config object for default 12C nucleus:")
pp(config_nucleus)
