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

#ifndef NNCrossSection_h
#define NNCrossSection_h
#include "nurex/Utils.h"
#include "nurex/Density.h"
#include "nurex/Nucleus.h"
#ifdef USE_THREADS
#include <mutex>
#include <atomic>
#endif

namespace nurex{

/// NNCrossSectionConstant
/// Returns constant nucleon-nucleon cross section, independent on energy
/// this us uusefull only for debug purposes
/// @param sigma_np set neutron-proton cross section
/// @param sigma_pp set proton-proton cross section
class NNCrossSectionConstant{
    public:
    NNCrossSectionConstant(double sigma_np, double sigma_pp):_np(sigma_np),_pp(sigma_pp){}
    double np(double)const{return _np;}
    double pp(double)const{return _pp;}

    private:
    double _np;
    double _pp;
};


///////////////////// NNCrossSectionFit /////////////////////////
///  NNCrossSectionFit
/// Fit to the measured proton-proton and neutron-proton cross-section data.
class NNCrossSectionFit{
    public:
    NNCrossSectionFit()noexcept{}
    NNCrossSectionFit(const NNCrossSectionFit&){}
    NNCrossSectionFit(NNCrossSectionFit&&){}
    ~NNCrossSectionFit(){}

    /// return np cross-section in mb
    double np(double energy);
    /// return nn and pp cross-section in mb
    double pp(double energy);

    double sigma_np(double energy) const;
    double sigma_pp(double energy) const;


    private:
    double energy_np=-1.0;
    double energy_pp=-1.0;
    double cs_np=0.0;
    double cs_pp=0.0;

    static double EXPO(double AA, double BB, double CC, double energy) noexcept;
    static double POLY(double AA, double BB, double CC, double DD, double energy) noexcept;
    static double POWER5(double AA, double BB, double CC, double DD, double EE, double FF, double energy) noexcept;
    #ifdef USE_THREADS
    std::mutex mutex;
    #endif
};

///////////////////// NNCrossSectionFile /////////////////////////
/// NNCrossSectionFile
/// Load proton-proton and neutron-proton cross-section from the file
/// It is assumed the files are txt format each with 2 columns.
/// 1st column corresponds to the enrgy in MeV/u, 2nd column to the measured cross-section in mb.
#ifndef NO_FILESYSTEM
class NNCrossSectionFile{
    public:

    /// Constructor
    /// @param np_filename - file with neutron-proton data points
    /// @param pp_filename - file with proton-proton data points
    NNCrossSectionFile(const char* np_filename="npcs.dat",const char* pp_filename="ppcs.dat");
    ~NNCrossSectionFile();
    NNCrossSectionFile(const NNCrossSectionFile& o){NNCrossSectionFile(o.npfile.c_str(), o.ppfile.c_str());}
    double np(double energy);
    double pp(double energy);

    private:
    std::string npfile;
    std::string ppfile;
    nurex::Interpolator *ipn;
    nurex::Interpolator *ipp;
};
#endif

/// FermiMotion
/// This add Fermi motion correction to the NNCrossSection
/// @tparam T - NNCrossSection class
///
/// Example (adding Fermi Motion correction to the NNCrossSectionFit):
/// \code{.cpp}
/// FermiMotion<NNCrossSectionFit> sigmaNN(90,90);
/// FermiMotion<NNCrossSectionFit> sigmaNN_halo(90,18);
/// \endcode
/// ~~~
template <typename T = NNCrossSectionFit>
class FermiMotion{
    public:
    /// Constructor
    /// @param p1 - Fermi Motion correction of the projectile in MeV/c
    /// @param p2 - Fermi Motion correction of the target in MeV/c
    FermiMotion(double p1=90.0, double p2=90.0);
    FermiMotion(const FermiMotion& o):FermiMotion(){p_variance = o.p_variance;p_sd = o.p_sd;}
    ~FermiMotion();

    /// change fermi motion corrections
    void SetMomentum(double p1, double p2){
        p_variance = p1*p1 + p2*p2;p_sd = sqrt(p_variance);prev_e_np=0.0; prev_e_pp=0.0;
        }

    double np(double energy);
    double pp(double energy);
    double variance()const{return p_variance;}
    private:

    double prev_e_np=-1.;
    double prev_np=0;
    double prev_e_pp=-1.;
    double prev_pp=0;
    T sigma_nn;
    double p_variance=0.0;
    double p_sd = 0.0;
    integrator_fm_type *integrator;
    #ifdef USE_THREADS
	std::mutex fm_mutex;
	#endif
};


/// FermiMotion Constructor
/// @param p1 - mean fermi momentum of the projectile nucleons
/// @param p2 - mean fermi momentum of the target nucleons
template <typename T>
FermiMotion<T>::FermiMotion(double p1, double p2){
    p_variance = p1*p1 + p2*p2;
    p_sd = sqrt(p_variance);
    integrator = new integrator_fm_type();
}

template <typename T>
FermiMotion<T>::~FermiMotion(){
    delete integrator;
}

/// neutron-proton cross section
/// @param energy energy in MeV/u
/// @return neutron-proton reaction cross section in mb
template <typename T>
double FermiMotion<T>::np(double energy){

    #ifdef USE_THREADS
    std::lock_guard<std::mutex> lock(fm_mutex);
    #endif

    if(energy==prev_e_np)return prev_np;

    double res;
    double p0 = p_from_T(energy);
    double pdev = p_sd;
    double range = fm_correction_range_factor*pdev;
    auto f = [&](double x){
            return sigma_nn.sigma_np(T_from_p(x))*gaussian(x,p0,p_variance);
            };
    auto f2 = [&](double x){
            return sigma_nn.sigma_np(T_from_p(x));
            };

    if(p0>range){
        res = integratorGH.integrate(f2,p0,pdev,true);
        //res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision);
        }
    else{
        res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision,4);
        //  double lo = p0-range;
        //  double hi = p0+range;
        //  res = integrator->integrate_intervals(f,
        //      {{lo,-lo},{-lo,p0},{p0,hi}},
        //      0,fm_integration_precision);
    }


    prev_np = res;
    prev_e_np = energy;
    return res;
}

/// proton-proton cross section
/// @param energy energy in MeV/u
/// @return proton-proton reaction cross section in mb
template <typename T>
double FermiMotion<T>::pp(double energy){

    #ifdef USE_THREADS
    std::lock_guard<std::mutex> lock(fm_mutex);
    #endif
    if(energy==prev_e_pp)return prev_pp;

    double res;
    double p0 = p_from_T(energy);
    double pdev = p_sd;
    double range = fm_correction_range_factor*pdev;
    auto f = [&](double x){
            return sigma_nn.sigma_pp(T_from_p(x))*gaussian(x,p0,p_variance);
            };

    auto f2 = [&](double x){
            return sigma_nn.sigma_pp(T_from_p(x));
            };

    if(p0>range){
        res = integratorGH.integrate(f2,p0,pdev, true);
        //res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision);
        }
    else{
        res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision,4);

        //  double lo = p0-range;
        //  double hi = p0+range;
        //  res = integrator->integrate_intervals(f,
        //      {{lo,-lo},{-lo,p0},{p0,hi}},
        //      0,fm_integration_precision);
        }

    prev_pp = res;
    prev_e_pp = energy;
    return res;
}


/// \class FermiMotion
/// This add Fermi motion correction to the NNCrossSection
/// @tparam T - SNNCrossSection class
/// @param p1 - mean fermi momentum of the projectile nucleons
/// @param p2 - mean fermi momentum of the target nucleons
///
/// Example (adding Fermi Motion correction to the NNCrossSectionFit):
/// ~~~ {.cpp}
/// FermiMotion<NNCrossSectionFit> sigmaNN(90,90);
/// FermiMotion<NNCrossSectionFit> sigmaNN_halo(90,18);
/// ~~~
template <typename T = NNCrossSectionFit>
class FermiMotionD{
    public:
    FermiMotionD();
    FermiMotionD(const FermiMotionD&):FermiMotionD(){}
    ~FermiMotionD();
    double np(double energy, double p1=0.0, double p2=0.0);
    double pp(double energy, double p1=0.0, double p2=0.0);
    //double np(double energy){return sigma_nn.np(energy);}
    //double pp(double energy){return sigma_nn.pp(energy);}
    private:
    T sigma_nn;
    integrator_fm_type *integrator;
};

/// FermiMotion Constructor
/// @param p1 - mean fermi momentum of the projectile nucleons
/// @param p2 - mean fermi momentum of the target nucleons
template <typename T>
FermiMotionD<T>::FermiMotionD(){
    integrator = new integrator_fm_type();
}

template <typename T>
FermiMotionD<T>::~FermiMotionD(){
    delete integrator;
}

/// neutron-proton cross section
/// @param energy energy in MeV/u
/// @return neutron-proton reaction cross section in mb
template <typename T>
double FermiMotionD<T>::np(double energy, double p1, double p2){
    double res;
    double p0 = p_from_T(energy);
    double pvar = p1*p1 + p2*p2;
    if(pvar<=0.0)return sigma_nn.np(energy);
    double pdev = sqrt(pvar);
    double range = fm_correction_range_factor*pdev;

    auto f = [&](double x){
            return sigma_nn.sigma_np(T_from_p(x))*gaussian(x,p0,pvar);
            };
    auto f2 = [&](double x){
            return sigma_nn.sigma_np(T_from_p(x));
            };

    if(p0>range){
        res = integratorGH.integrate(f2,p0,pdev,true);
        //res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision);
        }
    else{
        //res = ig.integrate(f2,p0,pdev,true);
        res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision);
    }
    return res;
}

/// proton-proton cross section
/// @param energy energy in MeV/u
/// @return proton-proton reaction cross section in mb
template <typename T>
double FermiMotionD<T>::pp(double energy, double p1, double p2){
    double res;
    double p0 = p_from_T(energy);
    double pvar = p1*p1 + p2*p2;
    if(pvar<=0.0)return sigma_nn.pp(energy);
    double pdev = sqrt(pvar);
    double range = fm_correction_range_factor*pdev;

    auto f = [&](double x){
            return sigma_nn.sigma_pp(T_from_p(x))*gaussian(x,p0,pvar);
            };

    auto f2 = [&](double x){
            return sigma_nn.sigma_pp(T_from_p(x));
            };

    if(p0>range){
        res = integratorGH.integrate(f2,p0,pdev, true);
        //res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision);
        }
    else{
        //res = ig.integrate(f2,p0,pdev,true);
        res = integrator->integrate(f,p0-range,p0+range,0,fm_integration_precision);
        }
    return res;
}

using NNCrossSection_FermiMotion = FermiMotion<NNCrossSectionFit>;
using NNCrossSection_FermiMotionD = FermiMotionD<NNCrossSectionFit>;
}
#endif
