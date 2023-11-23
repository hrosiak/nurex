import matplotlib.pyplot as plt
import numpy as np
import pynurex as nurex

p1 = nurex.make_evaporation_parameters({"level_density":"GC_RIPL"})

# calculating evaporation of excited 11C nucleus, ie 11C after 1n removal from 12C
# printing for one excitation energy
energy = 20.0
r = nurex.evaporation_ratios(A=11,Z=6,Ex=energy,evaporation_parameters=p1)
print("Evaporation channels ratios at excitation energy = ", energy)
print("R_neutron = ", r.n.G)
print("R_proton = ", r.p.G)
print("R_2h = ", r.d.G)
print("R_3h = ", r.t.G)
print("R_3he = ", r.he3.G)
print("R_alpha = ", r.a.G)
print("R_imf = ", r.imf.G)

#lets plot for whole excitation energy ranges
energies = np.linspace(0,60,50) # from 0 to 60 MeV
Rn = np.zeros_like(energies)
Rp = np.zeros_like(energies)
Rd = np.zeros_like(energies)
Rt = np.zeros_like(energies)
Rhe3 = np.zeros_like(energies)
Ra = np.zeros_like(energies)
Rimf = np.zeros_like(energies)
for i,energy in enumerate(energies):
    r = nurex.evaporation_ratios(A=11,Z=6,Ex=energy,evaporation_parameters=p1)
    Rn[i] = r.n.G
    Rp[i] = r.p.G
    Rd[i] = r.d.G
    Rt[i] = r.t.G
    Rhe3[i] = r.he3.G
    Ra[i] = r.a.G
    Rimf[i] = r.imf.G

plt.figure()
plt.plot(energies, Rn, label="n")
plt.plot(energies, Rp, label="p")
plt.plot(energies, Rd, label="d")
plt.plot(energies, Rt, label="t")
plt.plot(energies, Rhe3, label="3He")
plt.plot(energies, Ra, label="4He")
plt.plot(energies, Rimf, label="IMF")
plt.ylabel("Channel ratio")
plt.xlabel("Excitation energy [MeV]")
plt.legend()
plt.grid()
plt.show()