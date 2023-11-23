import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np

# the followin are 2 equivalent ways to create nuclear density.
# first via config dictionary
df1 = nurex.make_density({"type":"fermi","parameters":[5,0.7,0.1]})
# directly using pynurex class
df2 = nurex.DensityFermi(r=5,z=0.7,w=0.1)

# to get density at certain radius
d0 = df1.density(0.0)
print("Density at center = ", d0)

# now calculate Rrms, both should be the same
r1 = nurex.r_rms(df1)
r2 = nurex.r_rms(df2)
print("Rrms = ",r1)
print("Rrms = ",r2)



# now plot density as a function of radius
r = np.arange(0,8,0.1)  # radius from 0 to 8
rho = [df1.density(_r) for _r in r]  # density as a function of r
plt.figure()
plt.plot(r,rho)
plt.grid()
plt.xlabel("radius [fm]")
plt.ylabel("density [fm^-3]")
plt.show()

