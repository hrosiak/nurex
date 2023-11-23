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

#ifndef GlauberModel_h
#define GlauberModel_h
#include "nurex/Nucleus.h"
#include "nurex/NNCrossSection.h"
#include "nurex/Utils.h"
#include "nurex/ModelUtils.h"
#include "nurex/Models.h"
#include <memory>
#include <algorithm>
#include <cmath>
#ifdef USE_THREADS
#include <thread>
#endif
#include <type_traits>
#include <utility>
#include <cstdint>

namespace nurex{

/*
 * @param gm - GlauberModel class
 * @param E - energy of the projectile
 * return total reaction cross section
 */
template<typename GM>
double SigmaR(GM &gm,double E);

/*
 * @param gm - GlauberModel class
 * @param E - energy of the projectile
 * return total charge changing cross section
 */
template<typename GM>
double SigmaCC(GM &gm,double E);


/*
 * @param gm - GlauberModel class
 * @param E - energy of the projectile
 * return total 1-neutron and multi-neutron removal cross section
 */
template<typename GM>
double SigmaXN(GM &gm,double E);

/*
 * @param gm - GlauberModel class
 * @param E - energy of the projectile
 * return total 1-neutron removal cross section
 */
template<typename GM>
double Sigma1N(GM &gm,double E);

/*
 * @param gm - GlauberModel class
 * @param E - energy of the projectile
 * return neutron removal cross sections as array
 */
template<typename GM>
removals_type SigmaINs(const GM &gm, double E);

template<typename GM>
EvaporationProbabilities n_removals_evaporation(const GM &gm, double E);

/**
 *
 * scaling factor for energy as in Phys. Rev. C82(2010)014609, Eq.16
 * @param E - energy in MeV/u
 * @return correction factor for \f$\sigma_CC\f$
 */
inline double sigma_cc_scaling_factor(double E){
    return 1.141-E*6.507e-5;
}

/**
 * Coulomb correction for impact parameter
 * \param @E - Energy in MeV/u
 * \param @cs - Reaction Cross-Section in mb
 */
double coulomb_correction_simple(const Nucleus &projectile, const Nucleus &target, double E, double cs);

/**
 * relativistic Coulomb correction for impact parameter
 * \param @E - Energy in MeV/u
 * \param @cs - Reaction Cross-Section in mb
 */
double coulomb_correction_relativistic(const Nucleus &projectile, const Nucleus &target, double E, double cs);

/**
 * @brief charge changing correction due to proton evaporation, mk1
 * @param Projectile
 * @param Target
 * @param snx - x-neutron removal cross section in mb
 * @return - charged particle evaporation cross section in mb
 */
double cc_evaporation_cor(Nucleus const &Projectile, const removals_type&, const EvaporationParameters&);

/**
 * @brief charge changing correction due to proton evaporation, mk1
 * @param Projectile
 * @param Target
 * @param snx - 1neutron removal cross section
 * @return - particle evaporation cross section
 */
double total_evaporation_cor(Nucleus const &Projectile, double snx, const EvaporationParameters&);


removals_type epax_xn_ratios(const Nucleus&  Projectile, const Nucleus& Target, const double norm=1.0);

/////////////// Glauber Model class ////////////////////////
/**
 * /brief GlaberModel class
 */
template<template<typename> typename gm,  typename sigmann_type=NNCrossSectionFit>
class GlauberModel: public gm<GlauberModel<gm, sigmann_type>> {
    public: 
    GlauberModel(Nucleus *proj, Nucleus *tar, double range=0.0):projectile(*proj),target(*tar),_range(range){prepare();}
    GlauberModel(Nucleus& proj, Nucleus& tar, range_t range):projectile(proj),target(tar){if(range==range_t::FiniteRange){_range.pp=0.39;_range.pn=0.39;}prepare();}
    GlauberModel(Nucleus& proj, Nucleus& tar, double range=0.0):projectile(proj),target(tar),_range(range){prepare();}
    GlauberModel(Nucleus& proj, Nucleus& tar, sigmann_type snn, double range=0.0):projectile(proj),target(tar),sigma_nn(snn),_range(range){prepare();}
    const Nucleus& Projectile()const {return projectile;}
    const Nucleus& Target() const {return target;}
    Nucleus ProjectileCopy()const {return projectile;}
    Nucleus TargetCopy() const {return target;}
    double SigmaR(double E){return nurex::SigmaR(*this,E);}
    double SigmaCC(double E){return nurex::SigmaCC(*this,E);}
    double SigmaXN(double E){return nurex::SigmaXN(*this,E);}
    double Sigma1N(double E){return nurex::Sigma1N(*this,E);}
    removals_type SigmaINs(double E){return nurex::SigmaINs(*this,E);}
    EvaporationProbabilities n_removals_evaporation(double E) {return nurex::n_removals_evaporation(*this, E);}

    /// Set Coulomb Correction type
    void SetCoulombCorrection(coulomb_correction_t type){coulomb_correction = type;}

    /// Set Charge Changing Correction type
    void SetCCCorrection(cc_correction_t type){cc_correction = type;}

    /// Set Evaporation Parameters
    void SetEvaporationParameters(EvaporationParameters par){evaporation_parameters = par;}
    
    /// Set Evaporation Parameters
    void SetExcitationFunctionType(excitation_function_type par){evaporation_parameters.excitation_function = par;}

    /// return nucleon-nucleon sigma class
    sigmann_type& get_sigma_nn(){return sigma_nn;}

    template<nucleon_t NNTYPE>
    double sigmann(double);

    /// return range parameter    
    double beta() const { return _range.pp;}
    
    double range(nucleon_t NNTYPE) const { 
        if (NNTYPE == nucleon_t::pp || NNTYPE == nucleon_t::nn) {
            return _range.pp;
        }
        else if (NNTYPE == nucleon_t::pn || NNTYPE == nucleon_t::np) {
            return _range.pn;
        }
        else{
            assert(false);
            return 0.0;
        }        
        }

    /// return range type
    range_t range_type() const {return (_range.is_zero())?range_t::ZeroRange:range_t::FiniteRange;}
    const range_parameters& get_range() const {return _range;}

    /// set the range parameter
    void SetRange(double rpp, double rpn=-1){
        if(rpn<0.0)rpn=rpp;
        if(rpp>=0.0){
            _range.pp = rpp;
            _range.pn = rpn;
            this->init();
            status = 0; // reset status to recalculate profiles
            }
        }
    void SetRange(range_parameters r){
        _range.pp=r.pp;
        _range.pn=r.pn;
    }
    void prepare();

    Nucleus projectile;
    Nucleus target;
    sigmann_type sigma_nn;
    z_integrated_type z_integrated;

    int8_t status=0;        ///< to keep status of calculation
    int8_t special = 0;     ///< this is to indicate there is a special density
    coulomb_correction_t coulomb_correction = coulomb_correction_t::none;       ///< coulomb correction type
    cc_correction_t cc_correction = cc_correction_t::none;                       ///< charge-changing correction type
    double Eprev = 0;     ///< previously calculated energy
    range_parameters _range;  ///< range parameter
    EvaporationParameters evaporation_parameters; ///< store of excitation parameters


    using model_type = gm<GlauberModel<gm, sigmann_type>>;
    using nncrosssection_type = sigmann_type;
    //static_assert(has_phase<gm>::value, "Model type is required for model config class");
    //static_assert(valid_model<amodel_type<gm>>::value,"1st template parameter error: class with defined profile function must be provided");
    static_assert(valid_NNCrossSection<sigmann_type>::value,
                "2nd template parameter error: NNCrossSection class is required, the class must provide: pp(double T) and np(double T) functions returning pp and np nucleon-nucleon cross sections");
};

template<template<typename> typename gm,  typename sigmann_type>
void GlauberModel<gm,sigmann_type>::prepare(){

    if(!(projectile) || !(target)){
        throw std::invalid_argument("GlauberModel class not properly initialized");
    }

    // check if there is Dirac Density involved
    special = special_t::NN;
    if(projectile.GetDensityProton().type() == density_type::dirac) special = special|special_t::pp_Dirac;
//    else if(projectile.Z() == 1)  special|special_t::pp_Zero;
    if(projectile.GetDensityNeutron().type() == density_type::dirac) special = special|special_t::pn_Dirac;
    if(target.GetDensityProton().type() == density_type::dirac) special = special|special_t::tp_Dirac;
    if(target.GetDensityNeutron().type() == density_type::dirac) special = special|special_t::tn_Dirac;
    z_integrated.calculate(projectile,target);
    this->init();
}

template<template<typename> typename gm,  typename sigmann_type>
template<nucleon_t NNTYPE>
double GlauberModel<gm,sigmann_type>::sigmann(double E){
    if constexpr (NNTYPE == nucleon_t::pp || NNTYPE == nucleon_t::nn) {
        return 0.1 * sigma_nn.pp(E);
    }
    else{
        return 0.1 * sigma_nn.np(E);
    }
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////
template<typename GM>
double T(GM& gm, double b, double E)
{
    return exp(-2.0*gm.X(b, E));
};

template<typename GM>
double Tcc(GM& gm, double b, double E)
{
    if constexpr(can_do_cc<GM>){
        double x = 0.0;
        x = gm.Xpp(b, E);
        x += gm.Xpn(b, E);
        return exp(-2.0*x);
    }
    else{
        return -1.0;
    }
}

template<typename GM>
double Txn(GM& gm, double b, double E)
{
    if constexpr(can_do_cc<GM>){
        double xp = 0.0;
        double xn = 0.0;
            xp = gm.Xpp(b, E);
            xp += gm.Xpn(b, E);
            xn = gm.Xnn(b, E);
            xn += gm.Xnp(b, E);
        return exp(-2.0*xp)*(1-exp(-2.0*xn));
    }
    else{
        return -1.0;
    }
}

template<typename GM>
double SigmaR(GM &gm,double E){
    double res = -1;
    if (gm.Projectile().A() == 1.0 && gm.Target().A() == 1.0) {
        NNCrossSectionFit sigmaNN;
        if (gm.Projectile().Z() == gm.Target().Z())return sigmaNN.pp(E);
        else return sigmaNN.np(E);
    }

    // check if complete recalculation is needed
    gm.Calculate(E);

    auto f = [&](double b) {
        double bc = b;
        if(gm.coulomb_correction==coulomb_correction_t::sommerfeld){
            double beta = beta_from_T(E);
            bc = b_coulomb(b, gm.Projectile().A(), gm.Projectile().Z(), gm.Target().Z(), beta);
        }
        return (1 - exp(-2.0*gm.X(bc, E))) * b;
    };

    double Rmax = std::max(gm.z_integrated.pp.max(), gm.z_integrated.pn.max()) + std::max(gm.z_integrated.tp.max(), gm.z_integrated.tn.max());
    assert(Rmax>0.0);
    double C = 10.0*2.0* nurex::PI;
    res = C * integrator_adaptive.integrate(f, 0, Rmax, b_integration_precision/C);    
    assert(res>=0.0);
    if(gm.coulomb_correction == coulomb_correction_t::classic)res *= coulomb_correction_simple(gm.Projectile(),gm.Target(),E,res);
    else if(gm.coulomb_correction == coulomb_correction_t::relativistic)res *= coulomb_correction_relativistic(gm.Projectile(),gm.Target(),E,res);        
    return res;
}

template<typename GM>
double SigmaCC(GM &gm, [[maybe_unused]] double E)
{
    if (gm.Projectile().A() == 1.0 && gm.Target().A() == 1.0) {
        NNCrossSectionFit sigmaNN;
        if (gm.Projectile().Z() == gm.Target().Z())return sigmaNN.pp(E);
        else return sigmaNN.np(E);
    }

    if constexpr(can_do_cc<GM>){
        double res = -1;
        // check if complete recalculation is needed
        gm.Calculate(E);

        auto f = [&](double b) {
            double bc = b;
            if(gm.coulomb_correction==coulomb_correction_t::sommerfeld){
                double beta = beta_from_T(E);
                bc = b_coulomb(b, gm.Projectile().A(), gm.Projectile().Z(), gm.Target().Z(), beta);
                }
            return (1 - exp(-2.0*( gm.Xpp(bc, E) + gm.Xpn(bc, E)))) * b;
            };

        double Rmax = std::max(gm.z_integrated.pp.max(), gm.z_integrated.pn.max()) + std::max(gm.z_integrated.tp.max(), gm.z_integrated.tn.max());
        assert(Rmax>0.0);
        double C = 10.0*2.0* nurex::PI;
        res = C * integrator_adaptive.integrate(f, 0, Rmax, b_integration_precision/C);

        if(gm.coulomb_correction == coulomb_correction_t::classic)res *= coulomb_correction_simple(gm.Projectile(),gm.Target(),E,res);
        else if(gm.coulomb_correction == coulomb_correction_t::relativistic)res *= coulomb_correction_relativistic(gm.Projectile(),gm.Target(),E,res);

        if(gm.cc_correction == cc_correction_t::PRC82)res *= sigma_cc_scaling_factor(E);
        if(gm.cc_correction == cc_correction_t::evaporation || gm.cc_correction == cc_correction_t::test){
            removals_type nx_cs;
            if(gm.evaporation_parameters.xn_cs == 0){
                nx_cs = SigmaINs(gm,E);
            }
            else{
                double snx = SigmaXN(gm,E);
                nx_cs = epax_xn_ratios(gm.Projectile(), gm.Target(), snx);
            }
            double cor = cc_evaporation_cor(gm.Projectile(), nx_cs, gm.evaporation_parameters);
            res += cor;
        }
        assert(!std::isnan(res));
        return res;
    }
    else{
        return -1.0;
    }
}

inline double SigmaR(GlauberModelType &gm,double E){
  return gm.SigmaR(E);
};

inline double SigmaCC(GlauberModelType &gm,double E){
  return gm.SigmaCC(E);
};

template<typename GM>
double SigmaXN(GM &gm, double E)
{
    if (gm.Projectile().A() == 1.0 && gm.Target().A() == 1.0) {
        NNCrossSectionFit sigmaNN;
        if (gm.Projectile().Z() == gm.Target().Z())return sigmaNN.pp(E);
        else return sigmaNN.np(E);
    }

    if constexpr(can_do_cc<GM>){
        double res=-1.0;
        // check if complete recalculation is needed
        gm.Calculate(E);

        auto f = [&](double b) {
            double bc = b;
            if(gm.coulomb_correction==coulomb_correction_t::sommerfeld){
                double beta = beta_from_T(E);
                bc = b_coulomb(b, gm.Projectile().A(), gm.Projectile().Z(), gm.Target().Z(), beta);
                }
            return (Txn(gm, bc, E)) * b;
            };

        double Rmax = std::max(gm.z_integrated.pp.max(), gm.z_integrated.pn.max()) + std::max(gm.z_integrated.tp.max(), gm.z_integrated.tn.max());
        assert(Rmax>0.0);
        double C = 10.0*2.0* nurex::PI;
        res = C * integrator_adaptive.integrate(f, 0, Rmax, b_integration_precision/C);

        if(gm.coulomb_correction == coulomb_correction_t::classic)res *= coulomb_correction_simple(gm.Projectile(),gm.Target(),E,res);
        else if(gm.coulomb_correction == coulomb_correction_t::relativistic)res *= coulomb_correction_relativistic(gm.Projectile(),gm.Target(),E,res);
        assert(!std::isnan(res));
        res *= gm.evaporation_parameters.n_removal_scaling;
        return res;

    }
    else {
        return -1.0;
    }
}

template<typename GM>
double SigmaIN(GM &gmo, double E, int n)
{
    if (gmo.Projectile().A() == 1.0 && gmo.Target().A() == 1.0) {
        NNCrossSectionFit sigmaNN;
        if (gmo.Projectile().Z() == gmo.Target().Z())return sigmaNN.pp(E);
        else return sigmaNN.np(E);
    }

    if constexpr(can_do_cc<GM>){
        if(n > gmo.Projectile().N() || n<=0)return 0.0;
        double res=-1.0;

        // now we clone model, but with neutron normalization 1.0
        Nucleus proj = gmo.ProjectileCopy();
        Nucleus target = gmo.TargetCopy();
        const int Np = gmo.Projectile().N();
        assert(Np>0);
        if(n>Np || n<=0){
            return 0.0;
        }
        proj.NormalizeNeutronDensity(1.0);
        GlauberModel<OLA, typename GM::nncrosssection_type> gm(proj, target); // now fast model to get ratios
        gm.SetCoulombCorrection(gmo.coulomb_correction);
        gm.Calculate(E);

        auto f = [&](double b) {
            double bc = b;
            if(gm.coulomb_correction==coulomb_correction_t::sommerfeld){
                double beta = beta_from_T(E);
                bc = b_coulomb(b, gmo.Projectile().A(), gmo.Projectile().Z(), gmo.Target().Z(), beta);
                }
            double xp = 0.0;
            double xn = 0.0;
            xp = gm.Xpp(bc, E);
            xp += gm.Xpn(bc, E);
            xn = gm.Xnn(bc, E);
            xn += gm.Xnp(bc, E);
            double P = (1-exp(-2.0*xn));
           return b*(exp(-2.0*xp)*pow(P,n)*pow(1-P,Np-n));
            };

        double Rmax = std::max(gm.z_integrated.pp.max(), gm.z_integrated.pn.max()) + std::max(gm.z_integrated.tp.max(), gm.z_integrated.tn.max());
        assert(Rmax>0.0);
        double C = 10.0*2.0* nurex::PI;
        res = C * integrator_adaptive.integrate(f, 0, Rmax, b_integration_precision/C);

        if(gm.coulomb_correction == coulomb_correction_t::classic)res *= coulomb_correction_simple(gmo.Projectile(),gmo.Target(),E,res);
        else if(gm.coulomb_correction == coulomb_correction_t::relativistic)res *= coulomb_correction_relativistic(gmo.Projectile(),gmo.Target(),E,res);
        assert(!std::isnan(res));
        res *= gmo.evaporation_parameters.n_removal_scaling;
        return Combination(Np,n)*res;
    }
    else {
        return -1.0;
    }
}

template<typename GM>
removals_type SigmaINs(const GM &gmo, double E)
{
    constexpr int NMAX = 6;
    std::array<double,NMAX> res;
    if (gmo.Projectile().A() == 1.0 && gmo.Target().A() == 1.0) {
        NNCrossSectionFit sigmaNN;
        if (gmo.Projectile().Z() == gmo.Target().Z())res[0] = sigmaNN.pp(E);
        else res[0] =  sigmaNN.np(E);
    }

    if constexpr(can_do_cc<GM>){
        // now we clone model, but with neutron normalization 1.0
        Nucleus proj = gmo.ProjectileCopy();
        Nucleus target = gmo.TargetCopy();
        const int Np = gmo.Projectile().N(); 
        assert(Np>0);       
        proj.NormalizeNeutronDensity(1.0);
        GlauberModel<OLA, typename GM::nncrosssection_type> gm(proj,  target); // now fast model to get ratios
        gm.SetCoulombCorrection(gmo.coulomb_correction);
        gm.Calculate(E);

        int n = 1; // this is for looping but needs to be referenced by lambda
        auto f = [&](double b) {
            double bc = b;
            if(gm.coulomb_correction==coulomb_correction_t::sommerfeld){
                double beta = beta_from_T(E);
                bc = b_coulomb(b, gmo.Projectile().A(), gmo.Projectile().Z(), gmo.Target().Z(), beta);
                }
            double xp = 0.0;
            double xn = 0.0;
            xp = gm.Xpp(bc, E);
            xp += gm.Xpn(bc, E);
            xn = gm.Xnn(bc, E);
            xn += gm.Xnp(bc, E);
            double P = (1-exp(-2.0*xn));
            return b*(exp(-2.0*xp)*pow(P,n)*pow(1-P,Np-n));
            };

        double Rmax = std::max(gm.z_integrated.pp.max(), gm.z_integrated.pn.max()) + std::max(gm.z_integrated.tp.max(), gm.z_integrated.tn.max());
        assert(Rmax>0.0);
        double C = 10.0*2.0* nurex::PI;

        for(n = 1; n<=NMAX; n++){
            if(n > Np || n<=0){
                res[n-1] = 0.0;
                continue;
                }            
            res[n-1] = Combination(Np,n)* C * integrator_adaptive.integrate(f, 0, Rmax, b_integration_precision/C);
            res[n-1] *= gmo.evaporation_parameters.n_removal_scaling;
            if(gm.coulomb_correction == coulomb_correction_t::classic)res[n-1] *= coulomb_correction_simple(gmo.Projectile(),gmo.Target(),E,res[n-1]);
            else if(gm.coulomb_correction == coulomb_correction_t::relativistic)res[n-1] *= coulomb_correction_relativistic(gmo.Projectile(),gmo.Target(),E,res[n-1]);
            assert(!std::isnan(res[n-1]));
        }
    }    
    return res;
}


template<typename GM>
double Sigma1N(GM &gm, [[maybe_unused]] double E)
{
    if (gm.Projectile().A() == 1.0 && gm.Target().A() == 1.0) {return 0.0;}

    if constexpr(can_do_cc<GM>){
        double res = -1;
        gm.Calculate(E);  // check if complete recalculation is needed
        auto nx_cs = SigmaINs(gm,E);        
        res = nx_cs[0];  // -1n removal
        if(gm.coulomb_correction == coulomb_correction_t::classic)res *= coulomb_correction_simple(gm.Projectile(),gm.Target(),E,res);
        else if(gm.coulomb_correction == coulomb_correction_t::relativistic)res *= coulomb_correction_relativistic(gm.Projectile(),gm.Target(),E,res);                
        if(gm.cc_correction == cc_correction_t::evaporation){
            double cor = total_evaporation_cor(gm.Projectile(), res, gm.evaporation_parameters);
            res -= cor;
            }
        return res;
        }
    else{
        return -1.0;
    }
}

template<typename GM>
EvaporationProbabilities n_removals_evaporation(const GM &gm, double E){
    EvaporationProbabilities res;           
    if(gm.cc_correction != cc_correction_t::evaporation) return res;
    auto par = gm.evaporation_parameters;
    double Ex = Emax(gm.Projectile(), par);    

    for(int i = 0; i<res.Ptot.size();i++){
        if( (gm.Projectile().N()>i+1)){
            double Pc = charge_evaporation_probability_total(gm.Projectile().A()-1-i, gm.Projectile().Z(), Ex, i+1, par);
            double P =  total_evaporation_probability(gm.Projectile().A()-1-i,gm.Projectile().Z(),Ex,i+1,par);
            auto all = evaporation_probabilities(gm.Projectile().A()-1-i,gm.Projectile().Z(),Ex,i+1,par);
            res.Ptot[i] = P;            
            res.Pch[i] = Pc;
            res.Pn[i] = all.n.G;
            res.Pp[i] = all.p.G;
            res.Pd[i] = all.d.G;
            res.Pt[i] = all.t.G;
            res.Phe3[i] = all.he3.G;
            res.Pa[i] = all.a.G;
            res.Pimf[i] = all.imf.G;
        }        
        else {
            res.Ptot[i] = -1.0;            
            res.Pch[i] = -1.0;
            res.Pn[i] = -1.0;
            res.Pp[i] = -1.0;
            res.Pd[i] = -1.0;
            res.Pt[i] = -1.0;
            res.Phe3[i] = -1.0;
            res.Pa[i] = -1.0;
            res.Pimf[i] = -1.0;
        }
    }
    return res;
}


/**
 * @brief calculates max excitation energy per hole parameter for GS excitation distribution function
 * 
 * @param Projectile 
 * @param par 
 * @return double max excitation energy per hole in MeV
 */
inline double Emax(const Nucleus &Projectile, const EvaporationParameters &par){
  double Ex = par.Emax;
  if(Ex==0. || Projectile.A()<=4)return 0.0;
  if(Ex<=1.0){
    double rho0 = Projectile.DensityNeutron(0.0);
    double C = (Ex>0.)?Ex:Emax_fe_ratio;
    Ex = C*fermi_energy(rho0, neutron_mass);
  }
  return Ex;
}

template<typename GM>
class GMCompound{
    public:
    GMCompound(Nucleus& proj, Compound& tar, range_t range=range_t::ZeroRange);
    int num() const {return gm.size();}
    const Nucleus& get_nucleus(int i)const {return target.get_nucleus(i);}
    double get_fraction(int i)const {return target.molar_fraction(i);}
    Nucleus projectile;
    Compound target;
    std::vector<GM> gm;
};

template<typename GM>
GMCompound<GM>::GMCompound(Nucleus& proj, Compound& tar,range_t range):
projectile(proj),target(tar)
{
    if(!projectile || target.ncomponents()==0){
        throw std::invalid_argument("GlauberModelCompund class not properly initialized");
    }

    for(auto& e : target.get_elements()){
        gm.push_back(GM(projectile, e.nucleus, range));
    };
}

template<typename GM>
double SigmaR(GMCompound<GM> &gm, double E){

    double sum = 0.0;
    double stn_sum = 0.0;
    for(int i=0;i<gm.num();i++){
        stn_sum += gm.get_fraction(i);
        sum += gm.get_fraction(i) * SigmaR(gm.gm[i],E);
    }
    return sum/stn_sum;
}

/////////////////////// model aliases ////////////////////////
using GlauberModelOLA = GlauberModel<OLA, NNCrossSectionFit>;
using GlauberModelOLA_FM = GlauberModel<OLA, NNCrossSection_FermiMotion>;

using GlauberModelOLA = GlauberModel<OLA, NNCrossSectionFit>;
using GlauberModelOLA_FM = GlauberModel<OLA,NNCrossSection_FermiMotion>;

using GlauberModelMOL = GlauberModel<MOL>;
using GlauberModelMOL_FM = GlauberModel<MOL,NNCrossSection_FermiMotion>;


using GlauberModelMOL4C = GlauberModel<MOL4C>;
using GlauberModelMOL4C_FM = GlauberModel<MOL4C,NNCrossSection_FermiMotion>;


//////////////////// glauber model makers ////////////////////////////////////
template<template<typename> typename T, typename SIGMANN=NNCrossSectionFit>
auto make_glauber_model(Nucleus& p, Nucleus& t, range_t range = range_t::ZeroRange){
    //if constexpr(std::is_pointer<T>::value){
        //return new GlauberModel<typename std::remove_pointer<T>::type, SIGMANN>(p,t, range);
        //}
    //else{
        return GlauberModel<T, SIGMANN>(p,t, range);
      //  }
}

template<template<typename> typename T, typename SIGMANN>
auto make_glauber_model(Nucleus& p, Nucleus &t, SIGMANN& snn, range_t range = range_t::ZeroRange){
    //if constexpr(std::is_pointer<T>::value){
//        return new GlauberModel<typename std::remove_pointer<T>::type, SIGMANN>(p,t,snn, range);
//        }
    //else {
    return GlauberModel<T, SIGMANN>(p,t, snn, range);
//}
}

} // namespace nurex
#endif
