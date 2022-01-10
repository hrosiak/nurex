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
gm = nurex.make_model(model_config)

# to calculate single energy reaction cross section:
rcs = nurex.sigma_cc(gm, 350.0)
print("Sigma_R at 350MeV/u is {} mb".format(rcs))

# we can calculate charge changing cross section for multiple energies
energies = np.arange(100,900,20)  # from 100 to 900 MeV/u in step of 20 MeV/u
cccs = nurex.sigma_cc(gm, energies)

print("Charge-changing cross sections: ")
for i in range(len(energies)):
  print("{} MeV/u : {} mb".format(energies[i], cccs[i]))

# and to plot it

plt.figure()
plt.plot(energies, cccs)
plt.xlabel("Energy [MeV/u]")
plt.ylabel("sigma_cc [mb]")
plt.grid()
plt.show()

