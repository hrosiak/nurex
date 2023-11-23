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

#ifndef Models_h
#define Models_h
#include "nurex/NNCrossSection.h"
#include "nurex/Utils.h"
#include "nurex/Phase_functions.h"
#include "nurex/GlauberModelBase.h"
#include <memory>
#include <algorithm>

#ifdef USE_THREADS
#include <thread>
#endif

namespace nurex{

/// Here Glauber Models using different approximations are defined


/**
 * Modified Optical Limit Approximation Finite Range Model
 */ 
template<typename GM>
struct MOL{
    phaseshift_mol_type phase;
    range_integrated_type range_integrated;

    void init(){
        auto& gm = static_cast<GM&>(*this);
        range_integrated.calculate(gm.projectile, gm.target, gm.beta());
    }

    struct MOL_xyintegral{
        template <nucleonmol_t NNTYPE>
        static double phase_function(GM &gm, double b, double E){
                return (gm.range_type() == range_t::FiniteRange)
                            ?mol_xyintegral<NNTYPE, range_t::FiniteRange>(gm,b,E)
                            :mol_xyintegral<NNTYPE, range_t::ZeroRange>(gm, b, E);
        }
    };

    void Calculate(double E){
        auto& gm = static_cast<GM&>(*this);
        if(gm.status && E==gm.Eprev)return; // skip if its already calculate or the same energy
        Calculate_mol<MOL_xyintegral>(gm, E);
        gm.status = 1;
        gm.Eprev = E;
    }

    double X(double b, double){
        return phase._X.eval(b);
    }
};

/**
 * Modified Optical Limit Approximation Zero Range Model
 */
template<typename GM>
struct MOL4C{
    phaseshift_4c_type phase;
    range_integrated_type range_integrated;

    void init(){
        auto& gm = static_cast<GM&>(*this);
        if(gm.get_range().is_same())
            range_integrated.calculate(gm.projectile, gm.target, gm.beta());
    }

    void Calculate(double E){
        auto& gm = static_cast<GM&>(*this); 
        if(gm.status && E==gm.Eprev)return; // skip if its already calculate or the same energy
        Calculate_4c3(gm, E);
        gm.status = 1;
        gm.Eprev = E;
    }

    template <nucleon_t NNTYPE>
    double phase_function(double b, double E){
        auto& gm = static_cast<GM&>(*this);   
        double sigmann = 0.5*gm.template sigmann<NNTYPE>(E);
        auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
        auto target_z = get_target<NNTYPE>(gm.z_integrated);    
        auto projectile_zr = get_projectile<NNTYPE>(gm.range_integrated);
        auto target_zr = get_target<NNTYPE>(gm.range_integrated);
        double beta = gm.range(NNTYPE);
        
        if(gm.range_type() == range_t::FiniteRange){
            return (gm.get_range().is_same())?
                    mol4c_xyintegral_constrange_NN(projectile_z, target_z,projectile_zr, target_zr,b,sigmann)
                    :mol4c_xyintegral_NN<range_t::FiniteRange>(projectile_z, target_z,b, sigmann, beta);
        }
        else{
            return mol4c_xyintegral_NN<range_t::ZeroRange>(projectile_z, target_z,b, sigmann, 0.0);
        }                    
    }

    template <nucleon_t NNTYPE>
    double phase_function_dirac(double b, double E){
        auto& gm = static_cast<GM&>(*this);   
        double sigmann = 0.5*gm.template sigmann<NNTYPE>(E); // to be consistent profile function normalisatoin to 0.5
        auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
        auto target_z = get_target<NNTYPE>(gm.z_integrated);
        auto rho = (is_ftype<DiracFunction>(projectile_z)) ? target_z : projectile_z; //assign normal density to rho
                return mol4c_xyintegral_dirac(rho,b,sigmann,gm.beta());
        }

    double X(double b, double){
        double x =phase.Xpp.eval(b);
        x += phase.Xpn.eval(b);
        x += phase.Xnp.eval(b);
        x += phase.Xnn.eval(b);
        return x;
    }

    double Xpp(double b, double){
        return phase.Xpp.eval(b);
    }

    double Xnn(double b, double){
        return phase.Xnn.eval(b);
    }

    double Xnp(double b, double){
        return phase.Xnp.eval(b);
    }

    double Xpn(double b, double){
        return phase.Xpn.eval(b);
    }

};

/**
 * Optical Limit Approximation Zero Range Model
 */ 
template<typename GM>
struct OLA{
    //#define USE_CRANGE
    phaseshift_4c_type phase;
    #ifdef USE_CRANGE
    range_integrated_type range_integrated;
    #endif

    void init(){
        #ifdef USE_CRANGE
        auto& gm = static_cast<GM&>(*this);
        range_integrated.calculate(gm.projectile, gm.target, gm.beta());
        #endif
    }

    void Calculate(double E){
        auto& gm = static_cast<GM&>(*this);
        if(gm.status)return; // skip if its already calculated, different energy is corrected    
        Calculate_4c3(gm, E);
        gm.status = 1;
        gm.Eprev = E;
    }

        template <nucleon_t NNTYPE>
        double phase_function(double b, double){
                auto& gm = static_cast<GM&>(*this);                
                auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
                auto target_z = get_target<NNTYPE>(gm.z_integrated);
                double beta = gm.range(NNTYPE);
                #ifndef USE_CRANGE
                return ola_xyintegral_NN(projectile_z, target_z, b, beta);                
                #else
                if(beta>0.0){                    
                    auto projectile_zr = get_projectile<NNTYPE>(gm.range_integrated);
                    auto target_zr = get_target<NNTYPE>(gm.range_integrated);
                    return ola_xyintegral_constrange_NN(projectile_z, target_z,projectile_zr, target_zr,b);
                    //return ola_xyintegral_NN(projectile_z, target_z, b, beta); 
                    //return (gm.get_range().is_same())?ola_xyintegral_constrange_NN(projectile_z, target_z,projectile_zr, target_zr,b):ola_xyintegral_NN(projectile_z, target_z, b, beta); 
                }
                else{
                    return ola_xyintegral_NN(projectile_z, target_z, b,beta);
                }    
                #endif                        
                }

        template <nucleon_t NNTYPE>
        double phase_function_dirac(double b, double){
            auto& gm = static_cast<GM&>(*this);    
            auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
            auto target_z = get_target<NNTYPE>(gm.z_integrated);
            auto rho = (is_ftype<DiracFunction>(projectile_z)) ? target_z : projectile_z; //assign normal density to rho
            return ola_xyintegral_dirac(rho,b, gm.range(NNTYPE));
        }

    double X(double b, double E){
        auto& gm = static_cast<GM&>(*this);
        double x = 0.1*0.5*gm.sigma_nn.pp(E)*gm.phase.Xpp.eval(b);
        x += 0.1*0.5*gm.sigma_nn.np(E)*gm.phase.Xpn.eval(b);
        x += 0.1*0.5*gm.sigma_nn.np(E)*gm.phase.Xnp.eval(b);
        x += 0.1*0.5*gm.sigma_nn.pp(E)*gm.phase.Xnn.eval(b);
        return x;
    }

    double Xpp(double b, double E){
        auto& gm = static_cast<GM&>(*this);
        return 0.1*0.5*gm.sigma_nn.pp(E)*gm.phase.Xpp.eval(b);
    }

    double Xnn(double b, double E){
        auto& gm = static_cast<GM&>(*this);
        return 0.1*0.5*gm.sigma_nn.pp(E)*gm.phase.Xnn.eval(b);
    }

    double Xnp(double b, double E){
        auto& gm = static_cast<GM&>(*this);
        return 0.1*0.5*gm.sigma_nn.np(E)*gm.phase.Xnp.eval(b);
    }

    double Xpn(double b, double E){
        auto& gm = static_cast<GM&>(*this);
        return 0.1*0.5*gm.sigma_nn.np(E)*gm.phase.Xpn.eval(b);
    }
};

template<typename GM>
struct OLA_FMD{
    phaseshift_4c_type phase;
    range_integrated_type range_integrated;
    fp_z_integrated_type fp_integrated;
    double fe_coefficient=0.55;
 
    void init(){
        auto& gm = static_cast<GM&>(*this);
        //range_integrated.calculate(gm.projectile, gm.target, gm.beta());
        fp_integrated.calculate(gm.projectile, gm.target);
    }

    template <nucleon_t NNTYPE>
    double phase_function(double b, double E){
            auto& gm = static_cast<GM&>(*this);
            return ola_fm_xyintegral<NNTYPE>(gm, b, E);
            /*
            return (gm.range_type() == range_t::FiniteRange)?
                    ola_fm_xyintegral_constrange<NNTYPE>(gm,b, E)
                    //ola_fm_xyintegral<NNTYPE>(gm, b, E)
                    :ola_fm_xyintegral<NNTYPE>(gm, b, E);*/
        }

    template <nucleon_t NNTYPE>
    double phase_function_dirac(double b, double E){
            auto& gm = static_cast<GM&>(*this);
            //auto projectile_z = get_projectile<NNTYPE>(gm.z_integrated);
            //auto target_z = get_target<NNTYPE>(gm.z_integrated);
            //auto rho = (is_ftype<DiracFunction>(projectile_z)) ? target_z : projectile_z; //assign normal density to rho
            return ola_fm_xyintegral_dirac<NNTYPE>(gm, b, E);
        }

    void Calculate(double E){
        auto& gm = static_cast<GM&>(*this);
        if(gm.status && E==gm.Eprev)return; // skip if its already calculated or the same energy
        Calculate_4c3(gm, E);
        gm.status = 1;
        gm.Eprev = E;
    }
    void SetFECoefficient(double v){fe_coefficient = v;}

    template<nucleon_t NNTYPE>
    double sigma_nn_density(double rp, double rt, double e){
        double snn=0.0;
        auto& gm = static_cast<GM&>(*this);
        auto p_fp = get_projectile<NNTYPE>(fp_integrated);
        auto t_fp = get_target<NNTYPE>(fp_integrated);
        double p1 = p_fp->eval(rp);
        double p2 = t_fp->eval(rt);
        if constexpr (NNTYPE == nucleon_t::pp || NNTYPE == nucleon_t::nn) {
            snn = 0.5*0.1 * gm.sigma_nn.pp(e, fe_coefficient*p1, fe_coefficient*p2);
        }
        else{
            snn = 0.5*0.1 * gm.sigma_nn.np(e, fe_coefficient*p1, fe_coefficient*p2);
        }

        return snn;
    }

    double X(double b, double){
        double x =phase.Xpp.eval(b);
        x += phase.Xpn.eval(b);
        x += phase.Xnp.eval(b);
        x += phase.Xnn.eval(b);
        return x;
    }

    double Xpp(double b, double){
        return phase.Xpp.eval(b);
    }

    double Xnn(double b, double){
        return phase.Xnn.eval(b);
    }

    double Xnp(double b, double){
        return phase.Xnp.eval(b);
    }

    double Xpn(double b, double){
        return phase.Xpn.eval(b);
    }
};
template<typename GM>
struct MOL4C_FMD{

    phaseshift_4c_type phase;
    range_integrated_type range_integrated;
//    fe_z_integrated_type fe_integrated;
    fp_z_integrated_type fp_integrated;
    double fe_coefficient=0.55;

    void init(){
        auto& gm = static_cast<GM&>(*this);
        //range_integrated.calculate(gm.projectile, gm.target, gm.beta());
        //fe_integrated.calculate(gm.projectile, gm.target);
        fp_integrated.calculate(gm.projectile, gm.target);
    }

    template <nucleon_t NNTYPE>
    double phase_function(double b, double E){
        auto& gm = static_cast<GM&>(*this);        
        return (gm.range_type() == range_t::FiniteRange)?
            //mol4_fm_xyintegral_constrange<NNTYPE>(gm,b, E)
            mol4_fm_xyintegral<NNTYPE, range_t::FiniteRange>(gm, b, E)
            :mol4_fm_xyintegral<NNTYPE, range_t::ZeroRange>(gm, b, E);
        }

    template <nucleon_t NNTYPE>
    double phase_function_dirac(double b, double E){
        auto& gm = static_cast<GM&>(*this);        
        return (gm.range_type() == range_t::FiniteRange)?mol4_fm_xyintegral<NNTYPE, range_t::FiniteRange>(gm, b, E):
                                                    mol4_fm_xyintegral<NNTYPE, range_t::ZeroRange>(gm, b, E);
        }

    void Calculate(double E){
        auto& gm = static_cast<GM&>(*this);
        if(gm.status && E==gm.Eprev)return; // skip if its already calculate or the same energy
        Calculate_4c3(gm, E);
        gm.status = 1;
        gm.Eprev = E;
    }

    void SetFECoefficient(double v){fe_coefficient = v;}

    double X(double b, double){
        double x =phase.Xpp.eval(b);
        x += phase.Xpn.eval(b);
        x += phase.Xnp.eval(b);
        x += phase.Xnn.eval(b);
        return x;
    }

    double Xpp(double b, double){
        return phase.Xpp.eval(b);
    }

    double Xnn(double b, double){
        return phase.Xnn.eval(b);
    }

    double Xnp(double b, double){
        return phase.Xnp.eval(b);
    }

    double Xpn(double b, double){
        return phase.Xpn.eval(b);
    }
};
template<typename GM>
struct MOL_FMD{
    phaseshift_mol_type phase;
    range_integrated_type range_integrated;
    fp_z_integrated_type fp_integrated;
    double fe_coefficient=0.55;

    void init(){
        auto& gm = static_cast<GM&>(*this);
        range_integrated.calculate(gm.projectile, gm.target, gm.beta());
        //fe_integrated.calculate(gm.projectile, gm.target);
        fp_integrated.calculate(gm.projectile, gm.target);
    }

    struct phase_detail{
            template <nucleonmol_t NNTYPE>
            static double phase_function(GM &gm, double b, double E){
                    return (gm.range_type() == range_t::FiniteRange)?mol_xyintegral_constrange_fm<NNTYPE, range_t::FiniteRange>(gm,b,E):mol_xyintegral_constrange_fm<NNTYPE, range_t::ZeroRange>(gm, b, E);
            }
        };

    void Calculate(double E){
        auto& gm = static_cast<GM&>(*this);       
        if(gm.status && E==gm.Eprev)return; // skip if its already calculate or the same energy
        Calculate_mol<phase_detail>(gm, E);
        gm.status = 1;
        gm.Eprev = E;
    }

    void SetFECoefficient(double v){fe_coefficient = v;}

    double X(double b, double){
        return phase._X.eval(b);
    }

};

} // namespace nurex
#endif


