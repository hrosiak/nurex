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
#ifndef Phase_functions_h
#define Phase_functions_h

#include "nurex/Utils.h"
#include "nurex/ModelUtils.h"
#ifdef USE_THREADS
#include <thread>
#endif
#include <algorithm>
#include <type_traits>
#include <tuple>

/// Here the functions to caclulate phase shift function are defined
namespace nurex{

#define IS_DIRAC(x,y) is_ftype<DiracFunction>((x)) || is_ftype<DiracFunction>((y))

////////////////// 4component mol /////////////////////////////////////
inline double mol4c_xyintegral_dirac(const Functional *rho, double b, double sigmann, double beta){
    if(beta>0.0){
        auto f2 = [&](double x, double y) {
            return rho->eval(distance(x, y, b, 0.0));
        };
        double integral_rho_t_finite_range = integratorGH2D(f2,0, beta, 0 ,beta)/(2.0*PI*beta*beta);
        return (1.0 - exp(-sigmann * integral_rho_t_finite_range));
    }
    else{
        return (1.0 * (1 - exp(-sigmann * rho->eval(b))));
    }
}

template<range_t RANGE>
inline double mol4c_xyintegral_NN(const Functional *projectile_z, const Functional *target_z, double b, double sigmann, [[maybe_unused]]double beta){

    auto fintegrand = [&](double _x, double _y) {
        if constexpr (RANGE == range_t::FiniteRange){
            auto f_profile_rho = [&](double x, double y) {
                return target_z->eval(distance(x, y, 0.0, 0.0));
            };
            double integral_rho_t_finite_range = 0.5*integratorGH2D(f_profile_rho,_x, beta, _y ,beta)/(PI*beta*beta);
            return projectile_z->eval(distance(_x, _y, b, 0.0)) * (1 - exp(-sigmann * integral_rho_t_finite_range));
        }
        else{
            return projectile_z->eval(distance(_x, _y, b, 0.0)) * (1 - exp(-sigmann * target_z->eval(distance(_x, _y, 0.0, 0.0))));
        }
    };

    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();

    double rmax = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);
    double res = 0.0;
    res = integrator2D(fintegrand, lo, b, 0.0, rmax);
    res += integrator2D(fintegrand, b, hi, 0.0, rmax);

    std::swap(projectile_z, target_z);

    lo = std::max(-rmaxp, b-rmaxt);
    hi = std::min(rmaxp, b+rmaxt);
    res += integrator2D(fintegrand, lo, b, 0.0, rmax);
    res += integrator2D(fintegrand, b, hi, 0.0, rmax);

    return res; // 0.5 * 2.0(from integration)
}

inline double mol4c_xyintegral_constrange_NN(const Functional* projectile_z,
                                     const Functional* target_z,
                                     const Functional* projectile_zr,
                                     const Functional* target_zr,
                                     double b,
                                     double sigmann
                                     ){
    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();

    auto fintegrand = [&](double _x, double _y) {
        return projectile_z->eval(distance(_x, _y, b, 0.0)) * (1 - exp(-sigmann * target_zr->eval(distance(_x, _y, 0.0, 0.0))));
    };

    double rmaxy = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);

    double res = 0.0;
    res = integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    std::swap(projectile_z, target_z);
    std::swap(projectile_zr, target_zr);

    lo = std::max(-rmaxp, b-rmaxt);
    hi = std::min(rmaxp, b+rmaxt);

    res += integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);
    
    return res; //0.5 * 2.0(from integration)
}

///////////////////////////////// OLA //////////////////////////////////////////////
inline double ola_xyintegral_dirac(const Functional *rho, double b, double beta){
    if(beta>0.0){
        auto f2 = [&](double x, double y) {
            return rho->eval(distance(x, y, b, 0.0));
        };
        return integratorGH2D(f2, 0.0 , beta, 0.0 , beta)/(2.0*PI*beta*beta);
    }
    else {
        return rho->eval(b);
    }
}

inline double ola_xyintegral_NN(const Functional *projectile_z, const Functional *target_z, double b, [[maybe_unused]]double beta){
    // nucleon-nucleon case
    auto fintegrand = [&](double _x, double _y) {

        if (beta>0.0){ //Finite Range
            auto f_profile_rho = [&](double x, double y) {
                return target_z->eval(distance(x, y, 0.0, 0.0));
            };
            double integral_rho_t_finite_range = integratorGH2D(f_profile_rho,_x, beta, _y ,beta)/(2.0*PI*beta*beta);
            return projectile_z->eval(distance(_x, _y, b, 0.0)) * integral_rho_t_finite_range;
        }
        else { // Zero Range
            return projectile_z->eval(distance(_x, _y, b, 0.0)) * target_z->eval(distance(_x, _y, 0.0, 0.0));
        }
    };
    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();
    double rmaxy = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);
    double res = 0.0;
    res = integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    return 2.0 * res; // 2.0(from integration)
};

inline double ola_xyintegral_constrange_NN(const Functional *projectile_z, const Functional *target_z, const Functional *projectile_zr, const Functional *target_zr, double b){
    // nucleon-nucleon case
    auto fintegrand = [&](double _x, double _y) {
            return projectile_z->eval(distance(_x, _y, b, 0.0)) * target_zr->eval(distance(_x,_y,0.,0.));
    };
    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();
    double rmaxy = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);
    double res = 0.0;
    res = integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    return 2.0 * res; // 2.0(from integration)
};

template<nucleon_t NNTYPE, class GM>
    double sigma_nn_density(GM &gm, double rp, double rt, double e){
        double snn=0.0;
//        auto& gm = static_cast<GM&>(*this);
        auto p_fp = get_projectile<NNTYPE>(gm.fp_integrated);
        auto t_fp = get_target<NNTYPE>(gm.fp_integrated);
        double p1 = p_fp->eval(rp);
        double p2 = t_fp->eval(rt);
        if constexpr (NNTYPE == nucleon_t::pp || NNTYPE == nucleon_t::nn) {
            snn = 0.5*0.1 * gm.sigma_nn.pp(e, gm.fe_coefficient*p1, gm.fe_coefficient*p2);
        }
        else{
            snn = 0.5*0.1 * gm.sigma_nn.np(e, gm.fe_coefficient*p1, gm.fe_coefficient*p2);
        }

        return snn;
    }

template <nucleon_t NNTYPE,range_t RANGE, class GM>
double mol4_fm_xyintegral(GM &gm, double b, double E){
    double res = 0;

    auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
    auto target_z = get_target<NNTYPE>(gm.z_integrated);

    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();
    double beta = gm.range(NNTYPE);
    if (IS_DIRAC(projectile_z,target_z)) {
        auto rho = (is_ftype<DiracFunction>(projectile_z)) ? target_z : projectile_z;
        if constexpr (RANGE == range_t::FiniteRange){
            auto f2 = [&](double x, double y) {
                double rt = distance(x, y, b, 0.0);
                return sigma_nn_density<NNTYPE>(gm, rt, rt, E)*rho->eval(rt);
            };
            double integral_rho_t_finite_range = 0.5*integratorGH2D(f2,0, beta, 0 ,beta)/PI/beta/beta;
            return (1 - exp(-integral_rho_t_finite_range));
        }
        else{
            return sigma_nn_density<NNTYPE>(gm, b,b, E)*rho->eval(b);
        }
    } // end of if DiracFunction

    // nucleon-nucleon case
    auto fintegrand = [&](double _x, double _y) {
           if constexpr (RANGE == range_t::FiniteRange){
               auto f_profile_rho = [&](double x, double y) {
                   return target_z->eval(distance(x, y, 0.0, 0.0));
               };
               double integral_rho_t_finite_range = 0.5*integratorGH2D(f_profile_rho,_x, beta, _y ,beta)/PI/beta/beta;
               double rp = distance(_x, _y, b, 0.0);
               double rt = distance(_x, _y, 0.0, 0.0);
               double snn = sigma_nn_density<NNTYPE>(gm, rp,rt, E);
               return projectile_z->eval(distance(_x, _y, b, 0.0)) * (1 - exp(-snn * integral_rho_t_finite_range));
           }
           else{
               double rp = distance(_x, _y, b, 0.0);
               double rt = distance(_x, _y, 0.0, 0.0);
               double snn = sigma_nn_density<NNTYPE>(gm, rp,rt, E);
               return projectile_z->eval(rp) * (1 - exp(-snn * target_z->eval(rt)));
           }
    };

    double rmaxy = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);

    res = integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    std::swap(projectile_z, target_z);

    lo = std::max(-rmaxp, b-rmaxt);
    hi = std::min(rmaxp, b+rmaxt);

    res += integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    res = 0.5 * res;
    return 2.0 * res; // 2.0(from integration)
}

template <nucleon_t NNTYPE,class GM>
double mol4_fm_xyintegral_constrange(GM &gm, double b, double E){
    double res = 0;
    
    auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
    auto target_z = get_target<NNTYPE>(gm.z_integrated);
    auto projectile_zr = get_projectile<NNTYPE>(gm.range_integrated);
    auto target_zr = get_target<NNTYPE>(gm.range_integrated);

    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();
    // nucleon-nucleon case
    auto fintegrand = [&](double _x, double _y) {
            double rp = distance(_x, _y, b, 0.0);
            double rt = distance(_x, _y, 0.0, 0.0);
            double snn = sigma_nn_density<NNTYPE>(gm, rp,rt, E);
            return projectile_z->eval(rp) * (1 - exp(-snn * target_zr->eval(rt)));
    };

    double rmaxy = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);

    res = integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    std::swap(projectile_z, target_z);
    std::swap(projectile_zr, target_zr);

    lo = std::max(-rmaxp, b-rmaxt);
    hi = std::min(rmaxp, b+rmaxt);

    res += integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    res = 0.5 * res;
    return 2.0 * res; // 2.0(from integration)
}

template <nucleon_t NNTYPE, class GM>
double ola_fm_xyintegral_dirac(GM &gm, double b, double E){
    auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
    auto target_z = get_target<NNTYPE>(gm.z_integrated);
    bool p_dirac  = (is_ftype<DiracFunction>(projectile_z))?true:false;
    auto rho = (p_dirac) ? target_z : projectile_z; //assign normal density to rho
    double beta = gm.range(NNTYPE);
    if(beta>0.0){ // Finite Range
    auto f_profile_rho = [&](double x, double y) {
            double r = distance(x, y, b, 0.0);
            double snn = gm. template sigma_nn_density<NNTYPE>(r,r, E);
            return snn*rho->eval(r);
        };
    return integratorGH2D(f_profile_rho, 0.0, beta, 0.0, beta)/(2.*PI*beta*beta);
    }
    else{ // Zero Range
        double snn = gm.template sigma_nn_density<NNTYPE>(b,b, E);
        return snn*rho->eval(b);
    }
}

template <nucleon_t NNTYPE, class GM>
double ola_fm_xyintegral(GM &gm, double b, double E){
    static_assert(valid_NNCrossSection_FM<decltype(GM::sigma_nn)>::value,"NNCrossSection class does not have non-constant fermi momentum support");
    double res = 0;

    auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
    auto target_z = get_target<NNTYPE>(gm.z_integrated);

    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();
    double beta = gm.range(NNTYPE);
    // nucleon-nucleon case
    auto fintegrand = [&](double _x, double _y) {
        double rp = distance(_x, _y, b, 0.0);
        double rt = distance(_x, _y, 0.0, 0.0);
        double snn = gm.template sigma_nn_density<NNTYPE>(rp,rt, E);
        if (beta>0.0){
            auto f_profile_rho = [&](double x, double y) {
                return target_z->eval(distance(x, y, 0.0, 0.0));
            };
            double integral_rho_t_finite_range = 0.5*integratorGH2D(f_profile_rho,_x, beta, _y ,beta)/(PI*beta*beta);
            return snn*projectile_z->eval(rp) * integral_rho_t_finite_range;
        }
        else { // Zero Range
            return snn*projectile_z->eval(rp) * target_z->eval(rt);
        }
    };

    double rmaxy = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);

    res = integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);
    
    return 2.0 * res; //2.0(from integration)
}

template <nucleon_t NNTYPE, class GM>
double ola_fm_xyintegral_constrange(GM &gm, double b, double E){
    //static_assert(has_sigma_nn<typename GM::model_type>::value,"for ola_fm profile function proper NNCrossSection class must be provided with FM capability");
    static_assert(valid_NNCrossSection_FM<decltype(GM::sigma_nn)>::value,"NNCrossSection class does not have non-constant fermi momentum support");
    double res = 0;

    auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
    auto target_z = get_target<NNTYPE>(gm.z_integrated);
    auto projectile_zr = get_projectile<NNTYPE>(gm.range_integrated);
    auto target_zr = get_target<NNTYPE>(gm.range_integrated);

    double rmaxp = projectile_z->max();
    double rmaxt = target_z->max();    
    // nucleon-nucleon case
    auto fintegrand = [&](double _x, double _y) {
            double rp = distance(_x, _y, b, 0.0);
            double rt = distance(_x, _y, 0.0, 0.0);
            double snn = gm. template sigma_nn_density<NNTYPE>(rp,rt, E);
            return snn*projectile_z->eval(rp) * target_zr->eval(rt);
    };

    double rmaxy = std::min(rmaxt, rmaxp);
    double lo = std::max(-rmaxt, b-rmaxp);
    double hi = std::min(rmaxt, b+rmaxp);

    res = integrator2D(fintegrand, lo, b, 0.0, rmaxy);
    res += integrator2D(fintegrand, b, hi, 0.0, rmaxy);

    return 2.0 * res; //2.0(from integration)
}

//////////////////// MOL classical average, decomposed /////////////////////////


template <nucleonmol_t NNTYPE, class GM>
inline const std::tuple<const Functional*, const Functional*, const Functional*> get_z_integrated_mol(const GM &gm){

   if constexpr(NNTYPE == nucleonmol_t::projectile_p) {
        return std::make_tuple(&gm.z_integrated.pp, &gm.z_integrated.tp, &gm.z_integrated.tn);
    }
    else if constexpr(NNTYPE == nucleonmol_t::projectile_n) {
        return std::make_tuple(&gm.z_integrated.pn, &gm.z_integrated.tp, &gm.z_integrated.tn);
    }
    else if constexpr(NNTYPE == nucleonmol_t::target_p) {
        return std::make_tuple(&gm.z_integrated.tp, &gm.z_integrated.pp, &gm.z_integrated.pn);
    }
    else if constexpr(NNTYPE == nucleonmol_t::target_n) {
        return std::make_tuple(&gm.z_integrated.tn, &gm.z_integrated.pp, &gm.z_integrated.pn);
    }
    else{
        assert(false);
        return std::make_tuple(nullptr, nullptr, nullptr);
    }
}

template <nucleonmol_t NNTYPE, class GM>
inline const std::tuple<const Functional*, const Functional*, const Functional*> get_fe_z_integrated_mol(const GM &gm){

   if constexpr(NNTYPE == nucleonmol_t::projectile_p) {
        return std::make_tuple(&gm.fe_z_integrated.pp, &gm.fe_z_integrated.tp, &gm.fe_z_integrated.tn);
    }
    else if constexpr(NNTYPE == nucleonmol_t::projectile_n) {
        return std::make_tuple(&gm.fe_z_integrated.pn, &gm.fe_z_integrated.tp, &gm.fe_z_integrated.tn);
    }
    else if constexpr(NNTYPE == nucleonmol_t::target_p) {
        return std::make_tuple(&gm.fe_z_integrated.tp, &gm.fe_z_integrated.pp, &gm.fe_z_integrated.pn);
    }
    else if constexpr(NNTYPE == nucleonmol_t::target_n) {
        return std::make_tuple(&gm.fe_z_integrated.tn, &gm.fe_z_integrated.pp, &gm.fe_z_integrated.pn);
    }
    else{
        assert(false);
        return std::make_tuple(nullptr, nullptr, nullptr);
    }
}

template <nucleonmol_t NNTYPE, class GM>
inline const std::tuple<double, double> get_sigma_nn_mol(GM &gm, double E){
    double snn1=0, snn2=0;
    if constexpr(NNTYPE==nucleonmol_t::projectile_p){
        snn1 = 0.1*0.5*gm.sigma_nn.pp(E);
        snn2 = 0.1*0.5*gm.sigma_nn.np(E);

    }
    else if constexpr(NNTYPE==nucleonmol_t::projectile_n){
        snn1 = 0.1*0.5*gm.sigma_nn.np(E);
        snn2 = 0.1*0.5*gm.sigma_nn.pp(E);
    }
    else if constexpr(NNTYPE==nucleonmol_t::target_p){

        snn1 = 0.1*0.5*gm.sigma_nn.pp(E);
        snn2 = 0.1*0.5*gm.sigma_nn.np(E);
    }
    else if constexpr(NNTYPE==nucleonmol_t::target_n){

        snn1 = 0.1*0.5*gm.sigma_nn.np(E);
        snn2 = 0.1*0.5*gm.sigma_nn.pp(E);
    }
    else{
        assert(false);
    }
    return std::make_tuple(snn1, snn2);
}

template <nucleonmol_t NNTYPE>
inline const std::tuple<double, double> get_range_mol(const range_parameters &range){
    double r1=0, r2=0;
    if constexpr(NNTYPE==nucleonmol_t::projectile_p){
        r1 = range.pp;
        r2 = range.pn;
    }
    else if constexpr(NNTYPE==nucleonmol_t::projectile_n){
        r1 = range.pn;
        r2 = range.pp;
    }
    else if constexpr(NNTYPE==nucleonmol_t::target_p){
        r1 = range.pp;
        r2 = range.pn;
    }
    else if constexpr(NNTYPE==nucleonmol_t::target_n){
        r1 = range.pn;
        r2 = range.pp;
    }
    else{
        assert(false);
    }
    return std::make_tuple(r1, r2);
}

inline double mol_xyintegral_dirac(const Functional *rho2p, const Functional *rho2n, double snn1, double snn2, double b, double beta1=0.0, double beta2=0.0){
    double e1, e2;
    if(is_ftype<DiracFunction>(*rho2p)){
        if(beta1>0.0)
            e1 = snn1*finite_range(b,beta1);
        else
            e1 = (b==0.0)?snn1:0.0;
    }
    else{
        e1 = (beta1>0.0)?snn1*rho2p->eval(b):snn1;
    }

    if(is_ftype<DiracFunction>(*rho2n)){
        if(beta2>0.0)
            e2 = snn2*finite_range(b,beta2);
        else
            e2 = (b==0.0)?snn2:0.0;
    }
    else{
        e2 = (beta2>0.0)?snn2*rho2n->eval(b):snn2;
    }
    return 0.5*(1-exp(-e1-e2));
}

template <nucleonmol_t NNTYPE, range_t RANGE, class GM>
double mol_xyintegral(GM &gm, double b, double E){
    double res = 0;    
    const auto[rho1, rho2p, rho2n] = get_z_integrated_mol<NNTYPE>(gm);
    const auto[snn1, snn2] = get_sigma_nn_mol<NNTYPE>(gm, E);
    const auto[beta1, beta2] = get_range_mol<NNTYPE>(gm.get_range());
    assert(rho1 != nullptr); assert(rho2p != nullptr);assert(rho2n != nullptr);
    
    // case when one of the nuclei is proton, ie DiracDistribution
    // then the Eq 19 from manual is used.    
    if (is_ftype<DiracFunction>(*rho1)) {
        return mol_xyintegral_dirac(rho2p, rho2n, snn1, snn2, b, beta1, beta2);
    } // end of if DiracFunction

    // nucleon-nucleon case
    auto fintegrand = [&](double _x, double _y) {
        double e1, e2;
        if constexpr(RANGE == range_t::FiniteRange){
            double r = distance(_x, _y, 0.0, 0.0);
            auto f_profile_rho2p = [&](double x, double y) {
                return rho2p->eval(distance(x, y, 0.0, 0.0));
            };

            auto f_profile_rho2n = [&](double x, double y) {
                return rho2n->eval(distance(x, y, 0.0, 0.0));
            };

            if(is_ftype<DiracFunction>(*rho2p)){
                e1 = snn1*finite_range(r,beta1);
            }
            else{
                //e1 = snn1*rho2p->eval(r);
                e1 = snn1*0.5*integratorGH2D(f_profile_rho2p,_x, beta1, _y ,beta1)/PI/beta1/beta1;
            }
            if(is_ftype<DiracFunction>(*rho2n)){
                e2 = snn2*finite_range(r,beta2);
            }
            else{
                //e2 = snn2*rho2n->eval(r);
                e2 = snn2*0.5*integratorGH2D(f_profile_rho2n,_x, beta2, _y ,beta2)/PI/beta2/beta2;
            }
            return  rho1->eval(distance(_x, _y, b, 0.0))*(1-exp(-e1-e2));
        }
        else{ // Zero Range
            e1 = snn1 * rho2p->eval(distance(_x, _y, 0.0, 0.0));
            e2 = snn1 * rho2n->eval(distance(_x, _y, 0.0, 0.0));
            return rho1->eval(distance(_x, _y, b, 0.0)) * (1 - exp(-e1-e2));
        }
    };

    double rmax_p = rho1->max();
    double rmax_t = std::max(rho2p->max(), rho2n->max());

    double rmax = std::min(rmax_t, rmax_p);
    
    double lo = std::max(-rmax_t, b-rmax_p);
    double hi = std::min(rmax_t, b+rmax_p);

    res = integrator2D(fintegrand, lo, b, 0.0, rmax);
    res += integrator2D(fintegrand, b, hi, 0.0, rmax);
    res = 0.5 * res;
    return 2.0 * res; // 2.0(from integration)
}



template <nucleonmol_t NNTYPE, range_t RANGE, class GM>
double mol_xyintegral_constrange(GM &gm, double b, double E){
    double res = 0;
    const Functional *rho1;
    const Functional *rho2p;
    const Functional *rho2n;

    std::tie(rho1, rho2p, rho2n) = get_z_integrated_mol<NNTYPE>(gm);
    const auto[snn1, snn2] = get_sigma_nn_mol<NNTYPE>(gm, E);
    assert(rho1 != nullptr); assert(rho2p != nullptr);assert(rho2n != nullptr);

    if constexpr(RANGE==range_t::FiniteRange){
        if constexpr(NNTYPE==nucleonmol_t::projectile_p || NNTYPE==nucleonmol_t::projectile_n){
            rho2p = &gm.range_integrated.tp;
            rho2n = &gm.range_integrated.tn;
            }
        else if constexpr(NNTYPE==nucleonmol_t::target_p || NNTYPE==nucleonmol_t::target_n){
            rho2p = &gm.range_integrated.pp;
            rho2n = &gm.range_integrated.pn;
            }
        }

    // case when one of the nuclei is proton, ie DiracDistribution
    // then the Eq 19 from manual is used.

    if (is_ftype<DiracFunction>(*rho1)) {
        return mol_xyintegral_dirac(rho2p, rho2n, snn1, snn2, b, gm.beta());
    } // end of if DiracFunction

    // nucleon-nucleon case

    auto fintegrand = [&](double _x, double _y) {
        double e1, e2;
        double r = distance(_x, _y, 0.0, 0.0);
        if(is_ftype<DiracFunction>(*rho2p)){
            e1 = snn1*finite_range(r,gm.beta());
        }
        else{
            e1 = snn1*rho2p->eval(r);
        }
        if(is_ftype<DiracFunction>(*rho2n)){
            e2 = snn2*finite_range(r,gm.beta());
        }
        else{
            e2 = snn2*rho2n->eval(r);
        }
        double sum = rho1->eval(distance(_x, _y, b, 0.0))*(1-exp(-e1-e2));

        return sum;
    };

    double rmax_p = rho1->max();
    double rmax_t = std::max(rho2p->max(), rho2n->max());

    double rmax = std::min(rmax_t, rmax_p);

    double lo = std::max(-rmax_t, b-rmax_p);
    double hi = std::min(rmax_t, b+rmax_p);

    res = integrator2D(fintegrand, lo, b, 0.0, rmax);
    res += integrator2D(fintegrand, b, hi, 0.0, rmax);
    res = 0.5 * res;
    return 2.0 * res; // 2.0(from integration)
}

template <nucleonmol_t NNTYPE, range_t RANGE, class GM>
double mol_xyintegral_constrange_fm(GM &gm, double b, double E){
    double res = 0;
    const Functional *rho1;
    const Functional *rho2p;
    const Functional *rho2n;
    const Functional *fe1;
    const Functional *fe2p;
    const Functional *fe2n;


    std::tie(rho1, rho2p, rho2n) = get_z_integrated_mol<NNTYPE>(gm);
    std::tie(fe1, fe2p, fe2n) = get_z_integrated_mol<NNTYPE>(gm);
    assert(rho1 != nullptr); assert(rho2p != nullptr);assert(rho2n != nullptr);
    assert(fe1 != nullptr); assert(fe2p != nullptr);assert(fe2n != nullptr);

    if constexpr(RANGE==range_t::FiniteRange){
        if constexpr(NNTYPE==nucleonmol_t::projectile_p || NNTYPE==nucleonmol_t::projectile_n){
            rho2p = &gm.range_integrated.tp;
            rho2n = &gm.range_integrated.tn;
            }
        else if constexpr(NNTYPE==nucleonmol_t::target_p || NNTYPE==nucleonmol_t::target_n){
            rho2p = &gm.range_integrated.pp;
            rho2n = &gm.range_integrated.pn;
            }
        }

    auto sigmann = [&](double rp, double rt) -> std::tuple<double, double>{
        double s1=0.0,s2=0.0;
        double e1 = fe1->eval(rp);
        double e2p = fe2p->eval(rt);
        double e2n = fe2n->eval(rt);
        double p1 = (e1>0.0)?p_from_T(e1):0.0;
        double p2p = (e2p>0.0)?p_from_T(e2p):0.0;
        double p2n = (e2n>0.0)?p_from_T(e2n):0.0;

        if constexpr (NNTYPE == nucleonmol_t::projectile_p || NNTYPE == nucleonmol_t::target_p) {
            s1 = 0.5*0.1 * gm.sigma_nn.pp(E, gm.fe_coefficient*p1, gm.fe_coefficient*p2p);
            s2 = 0.5*0.1 * gm.sigma_nn.np(E, gm.fe_coefficient*p1, gm.fe_coefficient*p2n);
            }
        else{
            s1 = 0.5*0.1 * gm.sigma_nn.np(E, gm.fe_coefficient*p1, gm.fe_coefficient*p2p);
            s2 = 0.5*0.1 * gm.sigma_nn.pp(E, gm.fe_coefficient*p1, gm.fe_coefficient*p2n);
        }
        return {s1,s2};
    };

    // case when one of the nuclei is proton, ie DiracDistribution
    // then the Eq 19 from manual is used.

    if (is_ftype<DiracFunction>(*rho1)) {
        //return mol_xyintegral_dirac(rho2p, rho2n, snn1, snn2, b, gm.beta());
        return 0.0;
    } // end of if DiracFunction

    // nucleon-nucleon case

    auto fintegrand = [&](double _x, double _y) {
        double e1, e2;
        double rt = distance(_x, _y, 0.0, 0.0);
        double rp = distance(_x, _y, b, 0.0);
        auto [snn1, snn2] = sigmann(rp, rt);
        if(is_ftype<DiracFunction>(*rho2p)){
            e1 = snn1*finite_range(rt,gm.beta());
        }
        else{
            e1 = snn1*rho2p->eval(rt);
        }
        if(is_ftype<DiracFunction>(*rho2n)){
            e2 = snn2*finite_range(rt,gm.beta());
        }
        else{
            e2 = snn2*rho2n->eval(rt);
        }
        double sum = rho1->eval(rp)*(1-exp(-e1-e2));

        return sum;
    };

    double rmax_p = rho1->max();
    double rmax_t = std::max(rho2p->max(), rho2n->max());

    double rmax = std::min(rmax_t, rmax_p);

    double lo = std::max(-rmax_t, b-rmax_p);
    double hi = std::min(rmax_t, b+rmax_p);

    res = integrator2D(fintegrand, lo, b, 0.0, rmax);
    res += integrator2D(fintegrand, b, hi, 0.0, rmax);
    res = 0.5 * res;
    return 2.0 * res; // 2.0(from integration)
}

} // namespace nurex

#endif
