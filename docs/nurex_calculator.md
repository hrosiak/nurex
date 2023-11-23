Nurex Caluclator
================
Nurex Caluclator is a command line application and interface to the nurex library.

Usage
-----
The application is executed from the command line:

```
    nurex_calculator [options]";
```

the available arguments are:
 * --model MODEL_NAME - ie OLAZR, MOLFR_FM
 * --range value, optional for custom range parameter
 * --target - target symbol, ie 12c
 * --projectile - projectile symbol, ie 27al
 * --energy - energy in Mev/u, for energy range type Emin,Emax,Estep
 * --file FILENAME - save results to file FILENAME
 * --config jsonfile - load config from JSON formatted file


For simple calculation with default densities the calculator can be configured directly from command line, ie:
```
    nurex_calculator --model MOLZR_FM --target 12c --projectile 12c --energy 400;
    nurex_calculator --model OLAFR --target 12c --projectile 12c --energy 100,1000,10
```


For more complex configuration a json format input file can be used
```
    nurex_calculator --config c.json
```

Config File Format
------------------
The file must be a valid JSON formatted file.
The json file should contain the following keys: "model", "target", "projectile", "energy".


#### model
the __model__ key is string which specify the Glauber Model type, it corresponds to the command line argument __--model__
For available models see bellow.
example:
"model":"MOLFR"

for finite range models, the range parameter can be set to constant value using     the
__range__ keywordm which should take numerical value:

```
    "range":0.11
```

For models with the Fermi motion correction the momentum correction of projectile and target nucleons can be changed by __fermi_correction__ keyword.
The value is expected to be an array of two float numbers corresponding to MeV/c momentum corrections. The default value is 90.
```
    "fermi_motion":[85.0,85.0]
```

The following Coulomb Corrections types can be set by using the __coulomb_correction__ keyword, which accept string argument.
The following options are recognized:
 * __"none"__ - no Coulomb Correction is applied, this is default
 * __"classic"__
 * __"sommerfeld"__
 example:
 ```
    "coulomb_correction":"simple"
 ```


The following Charge-Changing Corrections types can be set by using the __charge_changing_correction__ keyword, which accept string argument.
The following options are recognized:
 * __"none"__ - no Coulomb Correction is applied, this is default
 * __"PRC82"__ - Energy dependent scaling factor as descrived in T. Yamaguchi, (2010)014609
 example:
 ```
    "coulomb_correction":"classic"
 ```



#### target and projectile
The __target__ and the __projectile__ keywords are either:
  * string - nucleus symbol, ie "12C", the default densities will be used
  * object - object with custom nucleus densities specified

It the target or projectile key is an obect it should contain the following keys: __"nucleus"__, __"proton_density"__ and __"neutron_density"__
  * __"nucleus"__ - string with nucleus symbol, ie "12C"
  * __"proton_density"__ - object with following format: {"type":"density_type","parameters":[parameter1, parameter2]}
  * __"neutron_density"__ - object with following format: {"type":"density_type","parameters":[parameter1, parameter2]}

the recognized densities with expected parameters format are listed below.

examples:
```
    "projectile":{
                    "nucleus":"64Ni"
                    "proton_density":{"type":"fermi","parameters":[1.4,1.2]},
                    "neutron_density":{"type":"fermi","parameters":[1.4,1.2]}
                  },
    "target":"12C"

```

#### Energy
The __energy__ keyword can be
1.a number specifying the kinetic energy:
```
    "energy":"500.0"
```

2. array of numbers for multiple energies:
```
    "energy":[100,200,500,1000]
```

3. Object specifying minimum energy, maximum energy and energy step, to calculate multiple energies:
```
    "energy":{
        "min": 100,
        "max": 1000,
        "step": 10
    }
```
instead of "step" key the "num" can be specified for integer number of steps between min and max energy.

#### File Output
To write results to the txt file __filke__ keyword should be used specifying a filename, ie
```
    "file":"output.txt"
```


Density Types
--------------
The following density types are recognized:
  * "ho" - Harmonic Oscillator type - parameters keyword is array of 2 numbers
  * "fermi" - 4 parameter fermi function - parameters keyword is array of 3 numbers, if array of 2 the 3rd parameter is set to 0.0
  * "dirac" - Dirac function - function is normalized
  * "file" - parameters keyword is string of filename with relative path to current directory, if file is not found parse error is returned

Models
-----------------
 * MOLFR - Modified Optical Limit
 * MOLFR_FM - Modified Optical Limit + Fermi Motion Correction
 * MOLZR - Modified Optical Limit + Zero Range Approximation
 * MOLZR_FM
 * OLAFR - Optical Limit Approximation
 * OLAFR_FM - Optical Limit Approximation + Fermi Motion Correction
 * OLAZR - Optical Limit Approximation + Zero Range Approximation
 * OLAZR_FM

Example File
-------------------
```
{
"model":"OLAZR_FM",
"target":"12C",
"projectile":"12C",
"energy":[350,450,600]
}

```



```
{
"model":"MOLFR",
"range":0.2,
"target":{
         "nucleus":"12c",
         "proton_density":{"type":"ho","parameters":[1.4,1.2]},
         "neutron_density":{"type":"ho","parameters":[1.4,1.2]}
        },

"projectile":{
            "nucleus":"64ni",
            "proton_density":{"type":"fermi_2p","parameters":[1.4,1.2]},
            "neutron_density":{"type":"fermi","parameters":[1.4,1.2]}
           },
"energy":{
    "min": 100,
    "max": 1000,
    "step": 10
    }
}

```
