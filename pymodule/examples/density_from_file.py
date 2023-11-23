import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

## Example how to load nuclear densities from file

########### 
# Method1:
#  Density can be loaded either by providing a table of density at radius
#  Array is assumed to be in format [[r_1, rho_1], [r_2, rho_2],...,[r_n, rho_n]]
#  Array can be file ie by loading from file

# parse file and save data to python list
rho_data = []
with open("48Ca.dat") as fr:
    for line in fr.readlines():
        r, density = line.split()
        rho_data.append( [float(r), float(density)] )

# making density from provided density list        
df = nurex.make_density({"type":"table","parameters":rho_data})

#########
# Method2:
# if file is formatted with 2 columns separated by whitespace filename can be directly provided
# Format of the file is assumed: 1st column is radius, 2nd column in density

df2 = nurex.make_density({"type":"File","parameters":"48Ca.dat"})


# now df and df2 are equivalent densities
# plotting densities

radii = np.arange(0,8,0.1)  # radius array from 0 to 8
rho = [df.density(r) for r in radii]  # densities to r
rho2 = [df2.density(r) for r in radii]  # densities to r

plt.figure()
plt.plot(radii,rho)
plt.plot(radii,rho2)
plt.grid()
plt.xlabel("radius [fm]")
plt.ylabel("density [fm^-3]")
#plt.show()


### Model Configuration with Densities from file

# using filename

model_config = {"model":"OLAZR_FM",
                  "target":"12C",
                  "projectile":{
                  "nucleus":"48Ca",
                  "proton_density":{"type":"file", "parameters":"48Ca.dat"},
                  "neutron_density":{"type":"file", "parameters":"48Ca.dat"}
                  },
                "charge_changing_correction":"evaporation",
                "level_density":"GC_GEM"
                  }


gm = nurex.make_model(model_config)

energies = np.arange(100,900,50)
cccs = nurex.sigma_cc(gm, energies)

plt.figure(figsize=(8,5))
plt.plot(energies, cccs)
plt.xlabel("Energy (MeV/u)")
plt.ylabel("charge-changing cross-section (mb)")
plt.grid()
plt.show()


model_config2 = {"model":"OLAZR_FM",
                  "target":"12C",
                  "projectile":{
                  "nucleus":"48Ca",
                  "proton_density":{"type":"table", "parameters":rho_data},
                  "neutron_density":{"type":"table", "parameters":rho_data}
                  },
                "charge_changing_correction":"evaporation",
                "level_density":"GC_GEM"
                  }

gm = nurex.make_model(model_config2)
cccs = nurex.sigma_cc(gm, 300)
print(cccs)
