/*
 *  Copyright(C) 2017, Andrej Prochazka, M. Takechi
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef NUCLEUS_H
#define NUCLEUS_H
#include "nurex/Density.h"
#include "nurex/NucID.h"
#include <vector>
#include <memory>
#include <utility>

namespace nurex{

/**
 * Nucleus definition, used for the arget and the projectile data
 *
 * example:
   \code{.cpp}
   Nucleus projectile1(62,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));
   auto density_proton = make_density<DensityFermi>(4.24,0.507);
   auto density_neutron = make_density<DensityFermi>(4.24,0.507);
   Nucleus projectile2(62,28,density_proton,density_neutron);
   Nucleus p_c("12c"); // get nucleus with default densities
   \endcode
 */
class Nucleus{
    public:
    Nucleus(){}

    /**
     * Constructor
     * @param A - mass number
     * @param Z - proton number
     * @param density_p - proton density
     * @param density_n - neutron density
     */
    Nucleus(int _A, int _Z, DensityType const &_density_p, DensityType const &_density_n);

    /**
     * Constructor - default Nucleus from the database
     * @param sym - nucleus symbol, ie "12C"
     */
    Nucleus(const char* sym);


    Nucleus(int _A, int _Z):a(_A),z(_Z){}
    Nucleus(const Nucleus &nuc);
    Nucleus(Nucleus &&nuc);

    Nucleus& operator=(Nucleus&& src);
    Nucleus& operator=(const Nucleus& )=delete;

    /// check if Nucleus is properly defined
    explicit operator bool()const;

    /// Reference to Proton Density Class
    const DensityType& GetDensityProton() const {return density_p;}

    /// Reference to Neutron Density Class
    const DensityType& GetDensityNeutron() const {return density_n;}

    /// density at radius r, density is sum of proton and neutron density
    /// @param r - radius in fm
    double Density(double r)const{return density_p.Density(r) + density_n.Density(r);}

    /// proton density at radius r
    /// @param r - radius in fm
    double DensityProton(double r)const{return density_p.Density(r);}

    /// neutron density at radius r
    /// @param r - radius in fm
    /// @return density
    double DensityNeutron(double r)const{return density_n.Density(r);}

    /// Rrms radii of the proton distribution
    /// @return Rrms in fm
    double RrmsProton()const {return density_p.Rrms();}

    /// Rrms radii of the neutron distribution
    /// @return Rrms in fm
    double RrmsNeutron()const {return density_n.Rrms();}

    void NormalizeProtonDensity(double n=1.0){density_p.Normalize(n);}
    void NormalizeNeutronDensity(double n=1.0){density_n.Normalize(n);}

    /// returns mass number
    int A() const {return a;}

    /// returns proton number
    int Z() const {return z;}

    /// returns neutron number
    int N() const {return a-z;}

    int a=-1;  /**< the mass number */
    int z=0;  /**< the proton number */

    /// check if 2 Nuclei are identical
    friend inline bool operator==(const Nucleus &n1, const Nucleus &n2);

    private:
    DensityType density_p;  /**< proton density */
    DensityType density_n;  /**< neutron density */
};
/// ////////////////

/**
 * get symbol of Nucleus
 * @param nuc - Nucleus reference
 * @return symbol or "-1" if Nucleus is not properly defined
 */
inline std::string symbol(const Nucleus &nuc){
    return std::to_string(nuc.A()) + get_element_symbol(nuc.Z());
}


///***** equality operator for Nucleus ***** ///
inline bool operator==(const Nucleus &n1, const Nucleus &n2){
    if(!n1 || !n2 || n1.A()!=n2.A() || n1.Z()!=n2.Z()){
        return false;
    }
    else{
        if( (n1.density_n != n2.density_n) || (n2.density_p != n2.density_p) ){
            return false;
        }
    }
    return true;
}

inline bool operator!=(const Nucleus &n1, const Nucleus &n2){
    return !(n1==n2);
}


struct atom{
        Nucleus nucleus;
        double stn;
    };


/**
 * Class used for compound targets
 */
class Compound{
    private:
    double molar_mass=0;
    std::vector<atom> atoms;
    public:
    Compound(std::initializer_list<atom>list);
    void add_element(const atom &a);
    unsigned int ncomponents() const {return atoms.size();}
    const Nucleus& get_nucleus(unsigned int i) const {return atoms[i].nucleus;}
    double molar_fraction(unsigned int i) const {return atoms[i].stn;}
    double M() const {return molar_mass;}
    std::vector<atom>& get_elements(){return atoms;}
};

} // nurex namespace
#endif
