#ifndef EVAPORATION_H
#define EVAPORATION_H
#include <cmath>
#include "nurex/numerics.h"
#include "nurex/data_frdm.h"
#include "nurex/ame16.h"
#include <optional>
#include <functional>
#include <stdexcept>
#include <map>

namespace nurex {

constexpr double Emax_fe_ratio = 1.0; // default

enum barrier_type {bass80=0,none=1, parametrized=2};
enum xn_type {gm=0, epax=1};
enum evaporation_preset_type {nurex=0, abla=1};
enum evaporation_config_type {
                              j_sum=1,
                              j_single = 2,
                              disable_k_factor = 4,
                              disable_imf=32,
                              disable_neutron=64,
                              test=128
                              };
enum level_density_type {
    GC_GEM = 0,
    GC_RIPL = 1,
    GC_KTUY05 = 2,
    ABLA = 255
};

enum excitation_function_type{
    GS = 0,
    CUSTOM = 1
};

struct EvaporationParameters{
    double Emax = -1.0;
    double n_removal_scaling = 1.0;
    uint8_t xn_cs = 0;
    uint8_t config = 0;
    uint8_t preset = 0;
    uint8_t density = level_density_type::GC_GEM;
    uint8_t excitation_function = excitation_function_type::GS;

#ifdef EVAPORATION_BARRIER_PARAMETRIZED
    int8_t barrier = barrier_type::parametrized;
#else
    int8_t barrier = 0;
#endif
};

struct emission_data{
    double G=0.;
    double rho=0.;
    double T=0.;
};

struct emission_results{
    emission_data g;
    emission_data n;
    emission_data p;
    emission_data d;
    emission_data t;
    emission_data he3;
    emission_data a;
    emission_data imf;
};

class ExcitationFunction{
    public:
    double w(double e, int a) const;
    void set(std::vector<double> energy, std::vector<double>rho, int a);
    void set_type(excitation_function_type t){type=t;}
    void reset(){ex.clear();norm.clear();}    
    double emax(int a)const {if(ex.find(a) != ex.end()){return ex.at(a).get_max();}else return 0.0;};
    std::map<int, Interpolator> ex;
    std::map<int, double> norm;
    excitation_function_type type;
};

extern EvaporationParameters default_evaporation;
extern ExcitationFunction custom_excitation;

double rho_ericson(double e, int a, const double Emax);
double w_ericson(double e, int a, const double Emax);
double cdf_w_ericson(double e, int a, const double Emax);

double rho_gs(double e, int a, const double Emax);
double w_gs(double e, int a, const double Emax);
double cdf_w_gs(double e, int a, const double Emax);
double cdf_wfx_gs(Functional& f, double e, int a, const double Emax);
double cdf_wfx_gs(const std::function<double(double)>&f, double e, int a, const double Emax, double hint=0.0);
double cdf_wfx_custom(const std::function<double(double)>&f, double e, int a);
double cdf_wfx(const std::function<double(double)>&f, double e, int a, const EvaporationParameters &config, const double Emax=0.0, double hint=0.0);


/**
 * @brief coulomb_potential of tow charged spheres
 * @param r - distance between centers
 * @param z1 - charge of 1st sphere
 * @param z2 - charge of 2nd sphere
 * @param r1 - radius of 1st sphere
 * @param r2 - radius of 2nd sphere
 * @return Coulomb potential in MeV
 */
double coulomb_potential(double r, int z1 , int z2, double r1, double r2);

/**
 * @brief bass - nuclear potential from bass80 model
 * @param r - distance between nuclei
 * @param A1 - 3rd square of atomic number of 1st nucleus
 * @param A2 - 3rd square of atomic number of 2nd nucleus
 * @return return nuclear potential in MeV
 */
double bass(double r, int A1, int A2);

/**
 * @brief nuclear_potential - this is sum of Coulomb + Bass potential
 * @param r
 * @param A1
 * @param Z1
 * @param A2
 * @param Z2
 * @return
 */
double nuclear_potential(double r, int A1, int Z1, int A2, int Z2);

/**
 * @brief fusion barrier - maximum of Coulomb - Bass barrier of 2 nuclei
 * @param A1
 * @param Z1
 * @param A2
 * @param Z2
 * @return barrier in MeV
 */
double fusion_barrier(int A1, int Z1, int A2, int Z2);

/**
 * @brief fusion barrier, parametrized from Kumari, CJP(2017)
 * @param A1
 * @param Z1
 * @param A2
 * @param Z2
 * @return barrier in MeV
 */
double fusion_barrier_parametrized(int A1, int Z1, int A2, int Z2);

/**
 * @brief charge_evaporation_probability
 * @param A - atomic number of prefragment nucleus
 * @param Z - proton number of prefragment nucleus
 * @param Emax - max Excitation energy per hole
 * @param h - number of holes
 * @return probability to evaporate charged particle
 */
double charge_evaporation_probability(int A, int Z, double Emax, int h, const EvaporationParameters &config = default_evaporation);
double charge_evaporation_probability_total(int A, int Z, double Emax, int h, const EvaporationParameters &config = default_evaporation);
double charge_evaporation_probability_simple(int A, int Z, double Emaxx, int h);

/**
 * @brief neutron_evaporation_probability
 * @param A - atomic number of prefragment nucleus
 * @param Z - proton number of prefragment nucleus
 * @param Emax - max Excitation energy per hole
 * @param h - number of holes
 * @return probability to evaporate neutron
 */
double neutron_evaporation_probability(int A, int Z, double Emax, int h, const EvaporationParameters &config = default_evaporation);


/**
 * @brief total particle evaporation probability
 * @param A - atomic number of prefragment nucleus
 * @param Z - proton number of prefragment nucleus
 * @param Emax - max Excitation energy per hole
 * @param h - number of holes
 * @return probability to evaporate particle
 */
double total_evaporation_probability(int A, int Z, double Emax, int h, const EvaporationParameters &config = default_evaporation);

/**
 * @brief total particle evaporation probability
 * @param A - atomic number of prefragment nucleus
 * @param Z - proton number of prefragment nucleus
 * @param Emax - max Excitation energy per hole
 * @param h - number of holes
 * @return probability to evaporate particle
 */
emission_results evaporation_probabilities(int A, int Z, double Ex, int a, const EvaporationParameters &config = default_evaporation);
/**
 * @brief T_freezeout
 * @param A - atomic number
 * @return freezeout temperature in MeV
 */
inline double T_freezeout(int A){
    return std::max(9.33*std::exp(-0.00282*A),5.5);
}

/**
 * @param A - mass number
 * @param a - level density parameter
 * @param Ex - effective excitation energy
 */
double shell_effect_damping(int A, double a, double Ex);

/**
 * @brief pairing energy
 * @param A
 * @param Z
 * @return
 */
double pairing_energy(int A, int Z, int conf=0);

/**
 * @brief a2_from_b2
 * @param b2 - beta_2 deformation parameter
 * @return
 */
inline double a2_from_b2(double b2){
    const double C = 5.0/(4.0*nurex::PI);
    return std::sqrt(C)*b2;
}

/**
 * @brief a4_from_b4
 * @param b4 - beta_4 deformation parameter
 * @return
 */
inline double a4_from_b4(double b4){
    const double C = 9.0/(4.0*nurex::PI);
    return std::sqrt(C)*b4;
}

/**
 * @brief bs_ratio - surface ratio deformed to spherical
 * @param b2 - b2 deforamtion parameter
 * @param b4 - b4 deforamtion parameter
 * @return ratio
 */
double bs_ratio(double b2, double b4);

/**
 * @brief bs_ratio - surface ratio deformed to spherical
 * @param A - mass number
 * @param Z - proton number
 * @return ratio
 */
inline double bs_ratio(int A, int Z){
    return bs_ratio(frdm::get_b2(A,Z), frdm::get_b4(A,Z));
}

/**
 * @brief asymptotic_level_density_parameter
 * @param A - mass of the    nucleus
 * @param bs - ratio between the surface deformed vs spherical nucleus
 * @return in MeV^-1
 */
double asymptotic_level_density_parameter(int A, double bs=1.0);

/**
 * @brief temperature_parameter
 * @param A - mass number
 * @param Z - proton number
 * @param Ex - Excitation energy, corrected effective
 * @param a - asymptotic level density
 * @return temperature
 */
inline double temperature_parameter(double E, double a){
    return std::sqrt(E/a);
};

/**
 * @param n - proton or neutron number
 * @return number of nucleons to the closest closed shell
 */
int closest_shell_difference(int n);

/**
 * @brief moments of inertia
 * @param A - mass number
 * @param b2 - b2 deformation parameter
 * @return pair of parallel and perpendicular moments of inertia
 */
std::pair<double, double> J(int A, double b2);

/**
 * @brief rotational enhancement
 * @param A - mass number
 * @param Z - proton number
 * @param Ex - excitation energy
 * @param sigma2_per - perpendicular spin cutoff
 * @param b2 - b2 deformation parameter
 */
double Krot(int A, int Z, double Ex, double sigma2_per, double b2);

/**
 * spin factor for level density
 */
double J_density_factor(double j, double sigma2);

struct prefragment{
    int A=0;
    int Z=0;
    double atilda=0.;
    double esp = 0.;
    double b2 = 0.;
    double pairing = -9999;

    double Sp = 0.;
    double Sn = 0.;
    double Sd = 0.;
    double St = 0.;
    double She3 = 0.;
    double Sa = 0.;

    double Cp = 0.0;
    double Cd = 0.;
    double Ct = 0.;
    double Che3 = 0.;
    double Ca = 0.;

    EvaporationParameters config;

    prefragment(){}
    prefragment(int A, int Z, const EvaporationParameters config=default_evaporation);
};

double S(prefragment &f, int Ap, int Zp);
double C(prefragment &f, int Ap, int Zp);

/**
 * decay width
 */
double width_e(prefragment &f, int Ap, int Zp, double B, double SB, double Ex, double j=0);
double width_gem(int Am, int Zm, int Ap, int Zp, double B, double SB, double Ex, double j=0);

double sigma_c(int Ad, int Ap, double epsilon, double Vb);

/**
 * @brief l_orb_distribution
 * @param Am - mass number of mother nucleus
 * @param Ad - mass number of daugter nucleus
 * @param l_m - orbital momentum of mother
 * @param Eeff - effective excitation energy
 * @param a - asymptotic level denisty
 * @return
 */
std::pair<double, double> l_orb_distribution(int Am, int Ad, double l_m, double Ef, double a);

/**
 * @brief angular_momentum_distribution
 * @param A - Mass number of projectile
 * @param Af - Mass number of fragment
 * @param beta - deformation parameter
 * @return pair of mean angular momentum of the fragment and sigma of the J distribution
 */
std::pair<double,double> angular_momentum_distribution(int A, int Af, double beta=0.0);

inline double mean_angular_momentum(int A, double beta=0.0){
    return sqrt(0.16*pow(A,2./3.)*(1.-(2.*beta/3.)));
}
double penetration_coefficient(int A, int Ap, double Temp);

/**
 * @brief Ecor - corrected energy used in Fermi-gas formula
 * @param f - prefragment structure
 * @param Ex - excitation energy to be corrected
 * @return
 */
double energy_corrected(const prefragment &f, double Ex);
double energy_corrected(int A, int Z, double Ex);

/**
 * @brief fermi_gas
 * @param A
 * @param Z
 * @param Ex - effective excitation energy, ie corrected for pairin
 * @return level density
 */
double fermi_gas_density(int A, int Z, double Ex);
double fermi_gas_density(const prefragment& f, double Ex);

/**
 * @brief const_temperature_density
 * @param f - prefragment
 * @param Ex - effective excitation energy
 * @return
 */
double const_temperature_density(const prefragment& f, double Ex);
/**
 * @brief temperature_parameter
 * @param f - prefragment
 * @param a - asymptotic level density
 * @return temperature
 */
double constant_temperature_parameter(const prefragment &f);

std::pair<double, double> level_density(prefragment &f, double Ex, double j=0);
inline std::pair<double, double> level_density(int A, int Z, double Ex, double j, const EvaporationParameters config=default_evaporation){
    prefragment f (A, Z);
    f.config = config;
    return  level_density(f, Ex, j);
};
double level_density_gem(prefragment &f, double E);
double level_density_ripl(prefragment &f, double E, double Em=-1);
double level_density_kawano(prefragment &f, double E, double Em=-1);

double superfluid_phase_critical_energy(const prefragment &f, double Ex);
double superfluid_phase_critical_energy(int A, int Z, double Ex);

emission_results evaporation_ratios(prefragment& f, double Ex, double j);
double charge_evaporation_function(prefragment& f, double Ex, double j=0.0);

double mean_decay_width(prefragment &d, int Ap, int Zp, double Ex, double lm);

emission_data get_emission_data(prefragment &f, int Ap, int Zp, double Ex, double jd = 0.0);
emission_data mean_emission_data(prefragment &f, int Ap, int Zp, double Ex, double jd);

// helpers
inline emission_results evaporation_ratios(int A, int Z, double Ex, double j=0, const EvaporationParameters config=default_evaporation){
    prefragment f(A,Z,config);
    return evaporation_ratios(f,Ex,j);
}

inline double decay_width(prefragment &f, int Ap, int Zp, double Ex, double jd){
    return get_emission_data(f,Ap, Zp, Ex, jd).G    ;
    }

inline double mean_decay_width(prefragment &f, int Ap, int Zp, double Ex, double lm){
    return mean_emission_data(f,Ap, Zp, Ex, lm).G;
    }

}// end of nurex namespace

#endif // EVAPORATION_H
