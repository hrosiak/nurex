import pynurex as nurex
import matplotlib.pyplot as plt
import numpy as np
import codecs


f=codecs.open ('projectile.txt', mode='r',encoding='utf-8')
line=f.readline()
list1=[]
while line:
      a=line.split()
      b=a[1]
      list1.append(b)
      line=f.readline()
f.close()

for projectilestr in list1:
       model_config = {"model":"OLAZR_FM", # model type
                  "projectile":projectilestr, # projectile, in this case default builtin 12C will be used
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

 # we can calculate charge changing cross section for multiple energies
       for j in range(100,200,50):
             energy = j
             rcs=gm.sigma_r(energy)
             cccs = gm.sigma_cc(energy)
             incs = gm.sigma_ins(energy)  # incs is array of primary neutron removals, incs[0] corresponds to 1n removal
             nevap = gm.n_removals_evaporation() #Total evaporation probability and charged-particle evaporation probabilities for i-neutron removals
 
             print(f"Energy = {energy} MeV/u")
             print(projectilestr)
             print(f"sigma_r = {rcs} mb,sigma_cc= {cccs} mb")
             path= 'results.txt'
             results=open(path,"a")
             results.write(f"projectile={projectilestr}  {energy} MeV/u  {rcs} mb  {cccs} mb"+'\n')
             results.close()
             for i in range(len(incs)):
                   Pt = nevap["Ptot"][i]
                   Pch = nevap["Pch"][i]
                   print(f"{i+1}n removal = {incs[i]} mb, Ptot = {Pt}, Pch = {Pch}")
                   results=open(path,"a")
                   results.write(f"{i+1}n {incs[i]} mb  {Pt}  {Pch}"+'\n')
                   results.close()

#f.close()


 
