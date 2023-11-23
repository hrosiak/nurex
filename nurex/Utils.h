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
#ifndef UTILS_H
#define UTILS_H

#include <functional>
#include <vector>
#include "nurex/Config.h"
#include "nurex/numerics.h"

namespace nurex{
    
template<typename T>
bool is_ftype(const Functional& density){
    return density.is_type<T>();
}

template<typename T>
bool is_ftype(const Functional* density){
    return density->is_type<T>();
}

/**
 * density distribution of finite range parameter
 * @param r - distance from center
 * @param b - range parameter in fm
 */
double finite_range(double r, double b);

/**
 * density distribution of finite range parameter, distance input is squared
 * @param r - distance squared from center, in fm^2
 * @param b - range parameter in fm
 */
double finite_range2(double r, double b);


/// calculates beta (velocity) from gamma factor
inline double gamma_from_beta(double beta){
    return 1.0/sqrt(1.0-beta*beta);
}

/// calculates gamma factor from beta (velocity)
inline double beta_from_gamma(double gamma){
    return sqrt(1.0-1.0/(gamma*gamma));
}

/// calculates momentum from velocity
inline double p_from_beta(double beta, double M=1.0){
    return M*atomic_mass_unit*beta*gamma_from_beta(beta);
}

/// calculates momentum from kinetic energy
inline double p_from_T(double T, double M=1.0){
    return M*sqrt(T*T + 2*T*atomic_mass_unit);
}

/// calculates kinetic energy from momentum
inline double T_from_p(double p, double M=1.0){
    return (sqrt((p*p) + (M*M*atomic_mass_unit*atomic_mass_unit))-(M*atomic_mass_unit))/M;
}

/// calculates gamma factor from kinetic energy
inline double gamma_from_T(double T){
    return (T+atomic_mass_unit)/atomic_mass_unit;
}

/// calculates beta from kinetic energy
inline double beta_from_T(double T){
    return beta_from_gamma(gamma_from_T(T));
}

/// calculates Ecm fom T, return in MeV units
inline double Ecm_from_T(double T, double Ap, double At){
    return T*Ap*At/(Ap+At);
}

/// calculates Ecm fom T, return in MeV units
inline double Ecm_from_T_relativistic(double T, double Ap, double At){
    double mp = Ap*atomic_mass_unit;
    double mt = At*atomic_mass_unit;
    double plab= p_from_T(T,Ap);
    double elab = sqrt(plab*plab + mp*mp);
    double ecm = sqrt(mp*mp + mt*mt + 2*elab*mt);
    //double pcm = plab * mt / ecm;
    //return sqrt(pcm*pcm+mp*mp)-mp;
    return elab + mt - ecm;
}

///Sommerfeld parameter
inline double sommerfeld_parameter(int z1, int z2, double beta){
return z1*z2*fine_structure_constant/beta;
}

inline double sommerfeld(int z1, int z2, double Ecm){
    return z1*z2*1.44*0.5/Ecm;
}

inline double closest_distance(double a1, int z1, int z2, double beta){
    return z1*z2*fine_structure_constant*hc_transp/(0.5*a1*atomic_mass_unit*beta*beta);
}

inline double b_coulomb(double b, double a1, int z1, int z2, double beta){
    double n = 0.5*closest_distance(a1, z1,z2,beta);
    return sqrt( (b*b) + (n*n) ) + n;
}

/**
 * @brief fermi energy
 * 
 * @param n density
 * @param mass nucleon mass in MeV/u
 * @return in MeV
 */
inline double fermi_energy(double n, double mass){
    return 0.5*hc_transp*hc_transp*pow(3*PI*PI*n,2.0/3.0)/mass;
}

/**
 * @brief fermi momentum
 * 
 * @param n density
 * @return momentum in MeV/c
 */
inline double fermi_momentum(double n){
    return hc_transp*pow(3*PI*PI*n,1.0/3.0);
}

inline double fermi_momentum_monitz(double A){
    return 259.416 - (152.842*std::exp(-A*9.5157e-2));
}

/*
 * Coulomb barrier in MeV, R1 and R2 are radii in fm unit
 */
inline double coulomb_barrier(int Z1, int Z2, double R1, double R2){
    return 1.44*Z1*Z2/(R1+R2);
}

/// gaussian function
inline double gaussian(double x, double mean,double variance){
    double dif = x-mean;
    return std::exp(-0.5*dif*dif/variance)/sqrt(2*PI*variance);
}


/// calculates the distance of the 2 points in 2D space
inline double distance(double x1, double y1, double x2, double y2)
{
    return sqrt(((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2)));
}

/// calculates the distance-squared of the 2 points in 2D space
inline double distance2(double x1, double y1, double x2, double y2)
{
    return ((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2));
}

/**
 * @brief makes vector of pair {J, P} J is angular momentum and P is probability
 * @param par - Gaussian parameters, pair of format {mean, width}
 * @return vector of pair of J and P distributed as gaussian
 */
std::vector<std::pair<int,double>> get_discrete_vector_by_area(std::pair<double,double> par);

/**
* @brief calculates number density per cm2 from provided thickness per cm2 and molar mass
*/
inline double number_density_cm2(double thickness_cm2, double molar_mass){
    return Avogadro*thickness_cm2/molar_mass;
    }

///////////////////////// type traits /////////////////////
template<class T, class=void>
struct has_sigma_nn : std::false_type{};
template<class T>
struct has_sigma_nn<T, std::void_t<decltype(T::sigma_nn)>> : std::true_type{};

template<class T, class=void>
struct valid_NNCrossSection : std::false_type{};
template<class T>
struct valid_NNCrossSection<T, std::void_t<decltype(std::declval<T>().np(0.0)),decltype(std::declval<T>().pp(0.0))>> : std::true_type{};

template<class T, class=void>
struct valid_NNCrossSection_FM : std::false_type{};
template<class T>
struct valid_NNCrossSection_FM<T, std::void_t<decltype(std::declval<T>().np(0.0,0.0,0.0)),decltype(std::declval<T>().pp(0.0,0.0,0.0))>> : std::true_type{};

template<class T, class=void>
struct has_4components : std::false_type{};
template<class T>
struct has_4components<T, std::void_t<decltype(std::declval<T>().Xpp(0.0, 0.0))>> : std::true_type{};
template<typename T>
constexpr auto can_do_cc = (has_4components<typename T::model_type>::value);

} // end of nurex namespace
#endif
