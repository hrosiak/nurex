{
"model":"MOLFR_FM",
"fermi_motion":[80,80],
"coulomb_correction":"relativistic",
"target":{
         "nucleus":"12c",
         "proton_density":{"type":"ho","parameters":[1.4,1.2]},
         "neutron_density":{"type":"ho","parameters":[1.4,1.2]}
        },
          
"projectile":{
            "nucleus":"64ni",
            "proton_density":{"type":"fermi","parameters":[3.5,0.7]},
            "neutron_density":{"type":"fermi","parameters":[3.5,0.7]}
           },
"energy":{
    "min": 100,
    "max": 1000,
    "step": 100
    }
}
