{
"model":"MOLFR_FM",
"fermi_motion":[80,80],
"target":{
         "nucleus":"12c",
         "proton_density":{"type":"ho","parameters":[1.4,1.2]},
         "neutron_density":{"type":"ho","parameters":[1.4,1.2]}
        },
          
    "projectile":{
            "nucleus":"48Ca",
            "proton_density":{"type":"file","parameters":"48Ca.dat"},
            "neutron_density":{"type":"file","parameters":"48Ca.dat"}
           },
"energy":{
    "min": 100,
    "max": 1000,
    "step": 100
    }
}
