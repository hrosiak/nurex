Config Format
-------------

Python dictionary or JSON format, key-value object,  is used to define the model configuration, projectile, target and densities.

## Model
keyword: __model__ - value is a string which specify the Glauber Model type. 
The model string constist of a combination of the followig 3 identifies:

Model types: 
 * OLA - Optical-Limit Approximation model
 * MOL - Modified Optical Limit model
 * MOL4C - 4 component MOL

Range types:
 * FR - Finite Range type
 * ZE - Zero Range approximation

Fermi motion correction

 * _FM - Fermi Motion correction 


Examples of valid __model__ strings:

 * MOLFR - Modified Optical Limit
 * MOLFR_FM - Modified Optical Limit + Fermi Motion Correction
 * MOLZR - Modified Optical Limit + Zero Range Approximation
 * MOLZR_FM
 * OLAFR - Optical Limit Approximation
 * OLAFR_FM - Optical Limit Approximation + Fermi Motion Correction
 * OLAZR - Optical Limit Approximation + Zero Range Approximation
 * OLAZR_FM

### Range
keyword: __range__ - numerical value of the range in fm for finite range model.
example:
```
    "range":0.39
```

### Fermi Motion Correction
keyword: __fermi_motion__ - array of 2 numerical values corresponding to the mean momentum correction in MeV/c. Defaults is [90,90]
example:
```
    "fermi_motion":[85.0,85.0]
```

### Fermi Energy Coefficient
keyword: __fermi_energy_coefficient__ - floating point value, default is 0.55, factor used to calculatate fermi motion correction from nucleon densities
example:
```
  "fermi_energy_coefficient":0.6
```

### Coulomb Correction
keyword: __coulomb_correction__ - string corresponding to the coulomb correction type. Allowed types are:
 * __"none"__ - no Coulomb Correction is applied, this is default
 * __"classic"__
 * __"sommerfeld"__ 

### Charge-Changing Correction
keyword: __charge_changing_correction__ - string corresponding to Charge-Changing Corrections type
The following options are recognized:
 * __"none"__ - no Coulomb Correction is applied, this is default
 * __"PRC82"__ - Energy dependent scaling factor as descrived in T. Yamaguchi, (2010)014609
 * __"evaporation"__
 
## target and projectile
The __target__ and the __projectile__ keywords are either:
  * string - nucleus symbol, ie "12C", the default densities will be used
  * object - object with custom nucleus densities specified

It the target or projectile key is an obect it should contain the following keys: __"nucleus"__, __"proton_density"__ and __"neutron_density"__
  * __"nucleus"__ - string with nucleus symbol, ie "12C"
  * __"proton_density"__ - object with following format: {"type":"density_type","parameters":[parameter1, parameter2]}
  * __"neutron_density"__ - object with following format: {"type":"density_type","parameters":[parameter1, parameter2]}

The following density types are recognized:
  * "ho" - Harmonic Oscillator type - parameters keyword is array of 2 numbers
  * "fermi" - 4 parameter fermi function - parameters keyword is array of 3 numbers, if array of 2 the 3rd parameter is set to 0.0
  * "dirac" - Dirac function - function is normalized
  * "gaussian" - Gaussian function - parameters is array of 1 number - width parameter
  * "table" - parameters is array of tabulated values of radius and density in format [[r0, rho0],[r2, rho1],...,[rn, rhon]] 

Example:
```
    "projectile":{
                    "nucleus":"64Ni"
                    "proton_density":{"type":"fermi","parameters":[1.4,1.2]},
                    "neutron_density":{"type":"fermi","parameters":[1.4,1.2]}
                  },
    "target":"12C"

``` 


Python example:
```python
model_config = {"model":"OLAZR_FM",
                  "projectile":{
                    "nucleus":"12c",
                    "proton_density":{"type":"ho", "parameters":[1.548415436,1.6038565]},
                    "neutron_density":{"type":"ho", "parameters":[1.548415436,1.6038565]}
                    },
                  "target":{
                  "nucleus":"12c",
                  "proton_density":{"type":"ho", "parameters":[1.548415436,1.6038565]},
                  "neutron_density":{"type":"ho", "parameters":[1.548415436,1.6038565]}
                  },
                "charge_changing_correction":"evaporation"
                  }
gm = make_model(model_config)
```

## Evaporation Parameters
### Emax
keyword: __excitation_energy__ - numerical value for the maximum excitation energy per hole in evaporation process. By default it is calculated from the provided neutron density.
Example:
```python
"excitation_energy": 40
```

### Level Densities
keyword: __level_density__ - string corresponding to the level density model
supported values:
  * "GC_GEM"
  * "GC_RIPL"
  * "GC_KTUY05"

### Neutron Removal Scaling
keyword: __n_removal_scaling__ - numerical value to scale calculated neutron removal cross sections. Can be usefull to adjust neutron removal cross section in evaporation corrections.
Example:
```python
"n_removal_scaling":0.75
```

### Disable Neutron Evaporation
keyword: __disable_neutron_evaporation__ - boolean, default False, to disable neutron evaporation channel

### Disable Intermediate Mass Fragments Evaporation
keyword: __disable_imf_evaporation__ - boolean, default False, to disable evaporation channels of particles with Z>2

