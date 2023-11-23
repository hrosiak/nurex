import sys
sys.path.insert(0,"../build")
import unittest
import pynurex as nurex
import math
import numpy as np

class TestStructures(unittest.TestCase):
    def test_make_density(self):
        df = nurex.DensityFermi(5,0.7)
        df2 = nurex.make_density({"type":"fermi","parameters":[5,0.7]})
        self.assertAlmostEqual(df.density(0.1),df2.density(0.1),4)
        self.assertAlmostEqual(df.density(3.1),df2.density(3.1),4)
        self.assertAlmostEqual(df.density(5.1),df2.density(5.1),4)

    def test_default_nucleus(self):
        nuc = nurex.default_nucleus("12c")
        nuc2 = nurex.default_nucleus(60,28)
        nuc3 = nurex.default_nucleus(6,1)

        self.assertTrue(nuc is not None)
        self.assertTrue(nuc2 is not None)
        self.assertTrue(nuc3 is None)

        self.assertTrue(nuc.A()==12)
        self.assertTrue(nuc.Z()==6)
        self.assertAlmostEqual(nuc.r_rms_proton(),2.332,1)
        self.assertAlmostEqual(nuc.r_rms_neutron(),2.332,1)

        self.assertTrue(nuc2.A()==60)
        self.assertTrue(nuc2.Z()==28)

    def test_density_type_string(self):
        df = nurex.DensityType(nurex.DensityFermi(5,0.7))
        nuc = nurex.default_nucleus("12c")
        self.assertEqual(nurex.density_type_to_string(df) , "fermi")
        self.assertEqual(nurex.density_type_to_string(nuc.get_density_proton()) , "ho")
        self.assertEqual(nurex.density_type_to_string(nuc.get_density_neutron()) , "ho")

    def test_nucleus_object(self):
        df1 = nurex.DensityType(nurex.DensityFermi(5,0.7))
        df2 = nurex.DensityType(nurex.DensityHO(1.5,1.2))
        nuc = nurex.Nucleus(13,6,df1,df2)
        r = nurex.nucleus_object(nuc)
        self.assertEqual(r["nucleus"], [13,6])
        self.assertEqual(r["proton_density"]["type"] , "fermi")
        self.assertEqual(r["proton_density"]["parameters"] , [5,0.7,0.0])
        self.assertEqual(r["neutron_density"]["type"] , "ho")
        self.assertEqual(r["neutron_density"]["parameters"] , [1.5,1.2])

    def test_model_object(self):
        config = {"model":"OLAFR_FM",
                  "projectile":{
                    "nucleus":"12c",
                    "proton_density":{"type":"ho", "parameters":[1.548415436,1.6038565]},
                    "neutron_density":{"type":"ho", "parameters":[1.548415436,1.6038565]}
                    },
                  "target":{
                  "nucleus":"12c",
                  "proton_density":{"type":"ho", "parameters":[1.548415436,1.6038565]},
                  "neutron_density":{"type":"ho", "parameters":[1.548415436,1.6038565]}
                  }
                  }
        gm = nurex.make_model(config)
        self.assertAlmostEqual(gm.sigma_r(30),1401,delta=5)
        self.assertAlmostEqual(gm.sigma_r(100),1010,delta=5)
        self.assertAlmostEqual(gm.sigma_r(200),837,delta=5)
        self.assertAlmostEqual(gm.sigma_cc(30),1265,delta=5)
        self.assertAlmostEqual(gm.sigma_cc(100),867,delta=5)
        self.assertAlmostEqual(gm.sigma_cc(200),689,delta=5)

        gm = nurex.make_model(config)
        self.assertAlmostEqual(gm.sigma_r(30),1401,delta=5)
        self.assertAlmostEqual(gm.sigma_r(100),1010,delta=5)
        self.assertAlmostEqual(gm.sigma_r(200),837,delta=5)

        energies = np.array([30,100,200])
        r = nurex.sigma_r(gm, energies)
        cc = nurex.sigma_cc(gm, energies)

        self.assertAlmostEqual(r[0],1401,delta=5)
        self.assertAlmostEqual(r[1],1010,delta=5)
        self.assertAlmostEqual(r[2],837,delta=5)
        self.assertAlmostEqual(cc[0],1265,delta=5)
        self.assertAlmostEqual(cc[1],867,delta=5)
        self.assertAlmostEqual(cc[2],689,delta=5)

    def test_range_param(self):
        config = {"model":"OLAFR",
                  "projectile":"12c",
                  "target":"12c"
                  }
        config2 = {"model":"OLA",
                  "projectile":"12c",
                  "target":"12c"
                  }
    
        config3 = {"model":"OLA", "target":"12c", "projectile":"12c", "range":[0.4, 0.5]}
        gm = nurex.make_model(config)
        gm2 = nurex.make_model(config2)
        gm3 = nurex.make_model(config3)
        e = 200
        self.assertTrue(gm.sigma_r(e)>gm2.sigma_r(e))
        self.assertEqual(gm2.get_range()["pp"], 0.0)
        self.assertEqual(gm2.get_range()["pn"], 0.0)
        gm2.set_range(0.39, 0.39)
        self.assertAlmostEqual(gm.sigma_r(e),gm2.sigma_r(e), delta=0.1)
        gm2.set_range(0.39,0.37)
        self.assertEqual(gm2.get_range()["pp"], 0.39)
        self.assertEqual(gm2.get_range()["pn"], 0.37)
        self.assertTrue(gm.sigma_r(e)>gm2.sigma_r(e))
        #self.assertTrue(gm.sigma_r(e)<gm3.sigma_r(e))
        gm2.set_range(0.4,0.5)
        self.assertEqual(gm2.get_range()["pp"], 0.4)
        self.assertEqual(gm2.get_range()["pn"], 0.5)
        self.assertEqual(gm3.get_range()["pp"], 0.4)
        self.assertEqual(gm3.get_range()["pn"], 0.5)
        self.assertTrue(gm.sigma_r(e)<gm2.sigma_r(e))
        self.assertAlmostEqual(gm2.sigma_r(e), gm3.sigma_r(e))



    def test_n_removal_scaling(self):
        config = {"model":"OLAFR_FM",
                  "projectile":"48Ca",
                  "target":"12C",
                  "charge_changing_correction":"evaporation",                  
                  }
        gm = nurex.make_model(config)
        config["n_removal_scaling"] = 0.7
        gm2 = nurex.make_model(config)
        config["charge_changing_correction"] = "none"
        config["n_removal_scaling"] = 1.0
        gm0 = nurex.make_model(config)

        dif = gm.sigma_cc(100) - gm0.sigma_cc(100)
        self.assertTrue(dif>0)
        self.assertTrue(gm.sigma_cc(100) > gm0.sigma_cc(100))        
        self.assertTrue(gm.sigma_cc(100) > gm2.sigma_cc(100))        
        self.assertTrue(gm2.sigma_cc(100) > gm0.sigma_cc(100))        
        self.assertAlmostEqual(gm2.sigma_cc(100)+dif*0.3, gm.sigma_cc(100), 1)

    def test_evaporation_parameters(self):
        conf = {"level_density":"GC_RIPL"}
        p1 = nurex.make_evaporation_parameters({})
        self.assertTrue(p1.density==0)
        p2 = nurex.make_evaporation_parameters(conf)
        self.assertTrue(p2.density==1)
        p3 = nurex.EvaporationParameters()
        self.assertTrue(p3.density==0)
        p3.density = nurex.level_density_type.GC_KTUY05
        self.assertTrue(p3.density==2)

    def test_evaporation_ratios(self):
        conf = {"level_density":"GC_RIPL"}
        p1 = nurex.make_evaporation_parameters(conf)
        r = nurex.evaporation_ratios(12,6,10,evaporation_parameters=p1)
        self.assertAlmostEqual(r.p.G,0.0,delta = 1e-3)
        self.assertAlmostEqual(r.g.G,0.0,delta = 1e-3)
        self.assertAlmostEqual(r.a.G,1.0,delta = 1e-3)
        conf = {"level_density":"GC_RIPL"}
        p1 = nurex.make_evaporation_parameters(conf)
        r = nurex.evaporation_ratios(A=12,Z=6,Ex=20,evaporation_parameters=p1)
        self.assertTrue(r.p.G>0.01)
        self.assertAlmostEqual(r.g.G,0.0,delta = 1e-3)
        self.assertTrue(r.a.G<0.9)
        p2 = nurex.make_evaporation_parameters({"level_density":"GC_GEM"})
        r2 = nurex.evaporation_ratios(A=12,Z=6,Ex=20,evaporation_parameters=p2)
        self.assertTrue(r2.p.G<r.p.G)
        self.assertTrue(r2.a.G>r.a.G)
        self.assertAlmostEqual(r2.he3.G,0.0,delta = 1e-3)
        r2 = nurex.evaporation_ratios(A=12,Z=6,Ex=38,evaporation_parameters=p2)
        self.assertTrue(r2.he3.G>0.0)
        self.assertTrue(r2.g.G==0.0)

    def test_neutron_removals_evaporation(self):
        model_config = {"model":"OLAZR_FM", # model type
                  "projectile":"12c", # projectile, in this case default builtin 12C will be used
                  "target":"12c",
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"classic", # or "none" or "sommerfeld", if not specified none is used
                "level_density":"GC_GEM"
                  }

        gm = nurex.make_model(model_config)
        nevap = gm.n_removals_evaporation() 
        for i in range(4):    
            pi = nevap["Ptot"][i]        
            self.assertTrue(pi<=nevap["Ptot"][i+1],f"i = {i}, Pi = {pi}")
            self.assertTrue(nevap["Pch"][i]<=nevap["Ptot"][i])
        self.assertAlmostEqual(nevap["Ptot"][2],1.0,delta = 1e-2)
        self.assertAlmostEqual(nevap["Pch"][2],1.0,delta = 1e-2)
    
    def test_custom_excitation(self):
        model_config = {"model":"OLAZR", # model type
                        "projectile":"12c", # projectile, in this case default builtin 12C will be used
                        "target":"12c",
                        "charge_changing_correction":"evaporation",
                        "coulomb_correction":"classic", # or "none" or "sommerfeld", if not specified none is used
                        "level_density":"GC_GEM"
                        }
        gm = nurex.make_model(model_config)
        model_config = {"model":"OLAZR", # model type
                        "projectile":"12c", # projectile, in this case default builtin 12C will be used
                        "target":"12c",
                        "charge_changing_correction":"none",
                        "coulomb_correction":"classic", # or "none" or "sommerfeld", if not specified none is used
                        "level_density":"GC_GEM"
                        }
        gm0 = nurex.make_model(model_config)


        par = nurex.EvaporationParameters()
        gm.set_excitation_function_type(nurex.excitation_function_type.GS)
        r1 = gm0.sigma_cc(1000)
        r2 = gm.sigma_cc(1000)
        
        ## by default custom excitation should be 0.0 so like no evaporation
        gm.set_excitation_function_type(nurex.excitation_function_type.CUSTOM)
        r3 = gm.sigma_cc(1000)
        self.assertTrue(r2>r3)
        self.assertAlmostEqual(r1, r3, delta=1)

        # now lets make custom excitation the same as GS
        energies = np.linspace(0,80,400)
        ex = nurex.Emax(gm.projectile(),gm.get_evaporation_parameters())
        self.assertAlmostEqual(ex,39.0,delta=5)        
        nurex.custom_excitation_reset()
        for i in range(1,5):   
            w = [nurex.w_gs(e,i,ex) for e in energies]
            nurex.custom_excitation_set(energies=energies, w=w, i=i)
            for j in range(100):
                self.assertAlmostEqual(nurex.custom_excitation_w(energy=energies[2*j],i=i), nurex.w_gs(energies[2*j],i,ex), delta=0.00001)    
        self.assertTrue(gm.get_evaporation_parameters().excitation_function==nurex.excitation_function_type.CUSTOM)
        r4 = gm.sigma_cc(1000)
        self.assertAlmostEqual(r2, r4, delta=2)

        # now squeezing GS curve to half range
        nurex.custom_excitation_reset()
        energies = np.linspace(0,40,400)
        for i in range(1,5):   
            w = [nurex.w_gs(2*e,i,ex) for e in energies]
            nurex.custom_excitation_set(energies=energies, w=w, i=i)
        r5 = gm.sigma_cc(1000.0)
        self.assertTrue(r5 < r2)
        self.assertTrue(r5 > r1)

if __name__ == "__main__":
    unittest.main()
