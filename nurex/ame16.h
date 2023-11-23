/*
 *  Copyright(C) 2017, Andrej Prochazka
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
#ifndef AME16_H
#define AME16_H

#include <cmath>

namespace ame16{

constexpr double electron_mass = 0.00054857991; //in amu unit
constexpr double amu = 931.494061; //in MeV
constexpr double proton_mass = 1.00727646692; // in amu
constexpr double neutron_mass = 1.0086649158; // in amu

/**
 * nucleus to id mapping
 */
constexpr inline int nucleus_id(int a, int z){
    return 10000*a + 10*z;
}

/**
 * @brief get ame 16 mass
 * @param nid - id of the nucleus
 * @return mass in amu unit
 */
double get_mass(int nid);

/**
 * @brief get ame 16 mass
 * @param a - mass number
 * @param z - proton number
 * @return mass in amu unit
 */
double inline get_mass(int a, int z){
    return get_mass(nucleus_id(a,z));
}

/**
 * Electron binding parametrization
 */
double inline Be(int z){
    return (14.4381*std::pow(z,2.39)) + (1.55468*1e-6*std::pow(z,5.35));
}

/**
 * @brief calculates nuclear mass
 * @param a - mass number
 * @param z - proton number
 * @return mass in amu unit
 */
double inline get_nuclear_mass(int a, int z){
    double mass = get_mass(nucleus_id(a,z));
    if(mass>0 && z>0){
        mass = mass - z*electron_mass + (Be(z)*1e-6/amu);
    }
    return mass;
}

/**
 * @brief Neutron(s) separation energy
 * @param a - mass number
 * @param z - proton number
 * @param n -number of removed neutrons, default is 1
 * @return n-neutron separation energy
 */
double inline Sn(int a, int z, int n=1){
    double s = 0.0;
    double mass = get_nuclear_mass(a,z);
    double mn = get_nuclear_mass(a-n,z);
    if(mass>0 && mn>0){
        s = -mass + mn + (n*neutron_mass);
        s *= amu;
    }
    return s;
}

/**
 * @brief Proton(s) separation energy
 * @param a - mass number
 * @param z - proton number
 * @param n -number of removed protons, default is 1
 * @return n-proton separation energy
 */
double inline Sp(int a, int z, int n=1){
    double s = 0.0;
    double mass = get_nuclear_mass(a,z);
    double mn = get_nuclear_mass(a-n,z-n);
    if(mass>0 && mn>0){
        s = -mass + mn + (n*proton_mass);
        s *= amu;
    }
    return s;
}

/**
 * @brief Alpha separation energy
 * @param a - mass number
 * @param z - proton number
 * @return alpha separation energy
 */
double inline Sa(int a, int z){
    double s = 0.0;
    double mass = get_nuclear_mass(a,z);
    double mn = get_nuclear_mass(a-4,z-2);
    double ma = get_nuclear_mass(4,2);
    if(mass>0 && mn>0){
        s = -mass + mn + ma;
        s *= amu;
    }
    return s;
}

/**
 * @brief separation energy to remove p number of protons and n number of neutrons
 * @param a - mass number
 * @param z - proton number
 * @param p - number of removed protons
 * @param n - number of removed neutrons
 * @return alpha separation energy
 */
double inline S(int a, int z, int p, int n){
    double s = 0.0;
    double mass = get_nuclear_mass(a,z);
    double mn = get_nuclear_mass(a-n-p,z-p);
    double ma = get_nuclear_mass(n+p,p);
    if(n==0)ma = p*proton_mass;
    if(p==0)ma = n*neutron_mass;
    if(mass>0 && mn>0 && ma>0){
        s = -mass + mn + ma;
        s *= amu;
    }
    return s;
}

/**
 * @brief Binding energy per nucleon
 * @param a - mass number
 * @param z - proton number
 * @return binding energy energy per nucleon
 */
double inline BA(int a, int z){
    double s = 0.0;
    double mass = get_nuclear_mass(a,z);
    if(mass>0 && a>1){
        s = -mass + (z*proton_mass) + ((a-z)*neutron_mass);
        s *= amu/a;
    }
    return s;
}


} //end of namespace
#endif
