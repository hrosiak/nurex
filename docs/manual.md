
Nucleon Density Distributions
------------------------------
In the NUREX, the following density types are defined and can be immediatelly used:
  * DensityFermi - 2- or 3-parameter Fermi type density
  * DensityHO - 2-oarameter Harmonic-Oscillator
  * DensityZero - Zero density, returns always 0
  * DensityDirac - Dirac type function
  * DensityGauss - Gaussian type function
  * DensityTable - tabulated density in form radius vs density

Example of defining density:
```cpp
DensityFermi d1(4.241,0.507);
DensityHO d2(1.548,1.603
``` 

The tabulated density can be specified using the std::vector<double> variables:
```cpp
std::vector<double> b_vector;
std::vector<double> density_vector;
... // load vectors here

// 3rd par. is optional normalization, defaults to 1.0
DensityTable dtable(b_vector, density_vector, 6.0); 
```
There is an helper function to load densities from a data file. The file is expected to have data in 2 columns, 1st column is radius parameter in fm unit, 2nd columns is the nucleon density in $fm^-3$ unit
```cpp
auto Ni64_p = density_from_file("Ni64p.dat",28)
```
### Nucleon Density 
The nucleon density at radius r can be calculated using __::Density(double r)__ function:
```cpp
    DensityFermi d1(4.241,0.507);
    double rho0 = d1.Density(0.0);
    double rho1 = d1.Density(1.0);
```

### Calculating root-mean-squared radius
The root-mean-squared radius of the distribution can be calculated as :
```cpp
    DensityFermi d1(4.241,0.507);
    double rrms = Rrms(d1);
```

Nucleus definition
-------------------
A Nucleus can be defined by providing A,Z and density distribution for protons and neutrons in form __Nucleus(A, Z, proton_density, neutron_density)
```cpp
Nucleus C12(12,6,DensityHO(1.54,1.6),DensityHO(1.55,1.6)); 

DensityFermi dp(4.24,0.51);
DensityFermi dn(4.26,0.5);
Nucleus Ni64(64,28, dp, dn);   
```

### Default nuclei
The nurex code provides default densities for many stable and unstable nuclei, for which charge radius is knows. The nucleus with default densities can be obtained as:
```cpp
Nucleus c12 = get_defaults_nucleus(12,6);
Nucleus c13 = get_defaults_nucleus("13c");
```

Glauber Model Types
--------------------
The model is created by specifying Glauber Model type and nucleon-nucleon cross section type. The Glauber Model class can be created specifying templated class __GlauberModel<MODEL, SIGMA\_NN>__  or  using the __make\_glauber\_model<MODEL, SIGMA\_NN>__ function:

```cpp
// 1. by specifying template parameters
// The 1st template argument is model type,
// the 2nd template argument is optional and 
// defines nn cross-section type, it defaults
// to NNCrossSectionFit class
// range is optional parameter, defaults to 0.0
GlauberModel<MODEL, SIGMA_NN> gm(projectile, target, range);

// 2. using make_glauber_model functions:
// for the template arguments the same apply as above
auto gm = make_glauber_model<MODEL, SIGMA_NN>(projectile, target)

// 3. if nn cross section class is supplied template parameter will be deduced:
auto gm = make_glauber_model<MODEL>(projectile, target, sigma_nn)
```
By default when SIGMA\_NN is not specified __NNCrossSectionFit__ class is used.

The following MODELs are available by default:
  * OLA
  * MOL
  * MOL4C

Also the following aliases can be used to directly select class type:
```cpp
using GlauberModelOLA = GlauberModel<OLA, NNCrossSectionFit>;
using GlauberModelOLA_FM = GlauberModel<OLA, NNCrossSection_FermiMotion>;

using GlauberModelOLA = GlauberModel<OLA, NNCrossSectionFit>;
using GlauberModelOLA_FM = GlauberModel<OLA,NNCrossSection_FermiMotion>;

using GlauberModelMOL = GlauberModel<MOL>;
using GlauberModelMOL_FM = GlauberModel<MOL,NNCrossSection_FermiMotion>;


using GlauberModelMOL4C = GlauberModel<MOL4C>;
using GlauberModelMOL4C_FM = GlauberModel<MOL4C,NNCrossSection_FermiMotion>;
```

### Selecting the range parameter
By default the model assume zero-range nn-interaction. This can be set using the __SetRange(double range)__ function. The range is in __fm__ unit.
```cpp
GlauberModelOLA gm(C12, C12);
gm.SetRange(0.39);
```

Examples of the model initialization:
```cpp
/// Optical-Limit Approximation with Zero-Range
GlauberModelOLA gm_zr(Ni64, C12);
GlauberModel<OLA> gm_zr2(Ni64,C12);

/// Optical-Limit Approximation with Finite-Range
// when range_t::FiniteRange, the range is 0.39
GlauberModel<OLA> gm_zr(Ni64,C12, range_t::FiniteRange);

//Modified Optical-Limit (MOL)
//Finite range
GlauberModel<MOL> gm_mol(Ni64,C12, range_t::FiniteRange);
GlauberModel<MOL> gm_mol1(Ni64,C12);
gm_mol1.SetRange(0.39);
```

