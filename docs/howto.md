How To
=======


Custom NN Cross Section
-----------------------

Custom NN Cross-section can be added by creating a class with functions: __np(double)__ and __pp(double)__. 
See sigmann_custom.cpp in examples directory.


*Example*: Custom NN Cross-section with returns factor of 0.9 of np() and 1.1 of pp() of supplied NNCrossSection, ie NNCrossSectionFit.

```cpp
class MyCrossSection{
    public:
    MyCrossSection(NNCrossSection *sigma):sigmann(sigma){};
    double np(double energy);
    double pp(double energy);

    NNCrossSection *sigmann;    
}

double MyCrossSection::np(double energy){
    double res;
    res = 0.9*sigmann->np(energy);
    return res;
}

double MyCrossSection::pp(double energy){
    double res;
    res = 1.1*sigmann->pp(energy);
    return res;
}

int main(){
    NNCrossSectionFit sigmann;
    MyCrossSection mysigma(&sigmann);
    ...
}
```


*Example:* NN Cross-section which returns constant value depending on the energy:

```cpp
class MyCrossSection{
    public:
    double np(double energy);
    double pp(double energy);
};
double MyCrossSection::np(double energy){
    if(energy>1000)
        return 40.0;
    else
        return 200.0;
}
double MyCrossSection::pp(double energy){
    if(energy>1000)
        return 45.0;
    else
        return 20.0;
}

int main(){
    MyCrossSection mysigma;
    ...
}
```
