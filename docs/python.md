Pynurex
=======
Nurex can be used in python using the C++ classes and functions bindigs and/or using a python dictionary configuration, where model configuration, projectile and target parameters are defined via python dictionary.

Python dictionary format 
========================
For specification how to define densities and model parametes see [JSON_format documentation](json_format.md).
The whole model can be created using the following function:
__pynurex.make_model(config:dict)__

Example:
```python
model_config = {"model":"OLAZR_FM",
                  "projectile":{
                    "nucleus":"12c",
                    "proton_density":{"type":"ho", "parameters":[1.548,1.603]},
                    "neutron_density":{"type":"ho", "parameters":[1.548,1.6038]}
                    },
                  "target":{
                  "nucleus":"12c",
                  "proton_density":{"type":"ho", "parameters":[1.548,1.6038]},
                  "neutron_density":{"type":"ho", "parameters":[1.548,1.6038]}
                  },
                "charge_changing_correction":"evaporation"
                  }
gm = make_model(model_config)
```

How To
==============

Change model range
------------------

Following will set range for pp part to 0.39fm and for pn to 0.37fm
```
config = {"model":"OLAFR",
           "projectile":"12c",
           "target":"12c"}
gm = nurex.make_model(config)
gm2.set_range(0.39,0.37)
```

