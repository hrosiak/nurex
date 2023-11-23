#ifndef ModelTypes_h
#define ModelTypes_h
#include "nurex/Nucleus.h"
#include "nurex/GlauberModelBase.h"
#include "nurex/Utils.h"
#include <optional>
#include <functional>
#include <algorithm>
#ifdef USE_THREADS
#include <thread>
#include <future>
#endif
namespace nurex{

/**
 *  returns z-integrated deinsity at impact parameter b
 */
double z_integral(const DensityType& density, double b);

/**
 * returns z-integrated densities cubic spline
 */
Functional ZIntegrate(const DensityType& nd);

/**
 * z- and range-integrated densities cubic spline
 * @param nd - Density
 * @param beta - range parameter
 */
Functional ZIntegrateRange(const DensityType& nd, double beta);


/**
 *  z-integrated fermi energy at impact parameter b
 * @param d - Density function
 * @param mass - mass of the nucleon
 * @param b - impact parameter
 */
double z_fermi_energy(const DensityType& d, double mass, double b);

/**
 *  z-integrated fermi momentum at impact parameter b
 * @param d - Density function
 * @param b - impact parameter
 */
double z_fermi_momentum(const DensityType& d, double b);

/**
 * z-integrated Fermi Energy cubic spline
 * @param nd - Density
 * @param mass - particle mass
 */
Functional ZIntegrate_Fermi_Energy(const DensityType& nd, double mass);

/**
 * z-integrated Fermi Momentum cubic spline
 * @param nd - Density
 * @param mass - particle mass
 */
Functional ZIntegrate_Fermi_Momentum(const DensityType& nd);

////////////////////////////////////////////////////////////
///////// z-integrated densities container //////////////////
////////////////////////////////////////////////////////////
struct z_integrated_type{
    Functional pp;
    Functional pn;
    Functional tp;
    Functional tn;
    void calculate(const Nucleus &projectile, const Nucleus &target){
        pp = (ZIntegrate(projectile.GetDensityProton()));
        pn = (ZIntegrate(projectile.GetDensityNeutron()));
        tp = (ZIntegrate(target.GetDensityProton()));
        tn = (ZIntegrate(target.GetDensityNeutron()));
    }
};

/////////////////////////////////////////////////////////////////
///////// Range integrated densities container //////////////////
/////////////////////////////////////////////////////////////////
struct range_integrated_type{
    Functional pp;
    Functional pn;
    Functional tp;
    Functional tn;
    double _beta=0.0;
    bool ready = false;
    void calculate(const Nucleus &projectile, const Nucleus &target, double beta){
        if( std::fabs(_beta-beta)<0.0001)return; // do not recalculate small difference in range
        if(beta==0.0){ // Zero Range case, make it empty
            pp = Functional();
            pn = Functional();
            tp = Functional();
            tn = Functional();
        }
        else{
            assert(beta>0.0);
            pp = (ZIntegrateRange(projectile.GetDensityProton(),beta));
            pn = (ZIntegrateRange(projectile.GetDensityNeutron(),beta));
            tp = (ZIntegrateRange(target.GetDensityProton(),beta));
            tn = (ZIntegrateRange(target.GetDensityNeutron(),beta));
        }
        _beta = beta;
        ready = true;
    }
};

////////////////////////////////////////////////////////////
///////// z-integrated Fermi Energies container ////////////
////////////////////////////////////////////////////////////
struct fe_z_integrated_type{
    Functional pp;
    Functional pn;
    Functional tp;
    Functional tn;
    void calculate(const Nucleus &projectile, const Nucleus &target){
        pp = (ZIntegrate_Fermi_Energy(projectile.GetDensityProton(), proton_mass));
        pn = (ZIntegrate_Fermi_Energy(projectile.GetDensityNeutron(), neutron_mass));
        tp = (ZIntegrate_Fermi_Energy(target.GetDensityProton(), proton_mass));
        tn = (ZIntegrate_Fermi_Energy(target.GetDensityNeutron(), neutron_mass));
    }
};

struct fp_z_integrated_type{
    Functional pp;
    Functional pn;
    Functional tp;
    Functional tn;
    void calculate(const Nucleus &projectile, const Nucleus &target){
        pp = (ZIntegrate_Fermi_Momentum(projectile.GetDensityProton()));
        pn = (ZIntegrate_Fermi_Momentum(projectile.GetDensityNeutron()));
        tp = (ZIntegrate_Fermi_Momentum(target.GetDensityProton()));
        tn = (ZIntegrate_Fermi_Momentum(target.GetDensityNeutron()));
    }
};

/*
 * return projectile part of T type for NNTYPE nucleon-nucleon combination
 * ie.
 * auto p_fe = get_projectile<nucleon_t::pp>(gm.fe_integrated);
 */
template <nucleon_t NNTYPE, typename T>
const Functional* get_projectile(const T &t){

    if constexpr (NNTYPE == nucleon_t::pp || NNTYPE ==  nucleon_t::pn) {
        return &t.pp;
    }
    else if(NNTYPE == nucleon_t::np || NNTYPE ==  nucleon_t::nn){
        return &t.pn;
    }    
    else{
        assert(false);
        return nullptr;
    }
}

template<typename T>
const Functional* get_projectile(const T &t, nucleon_t NNTYPE){
    if (NNTYPE == nucleon_t::pp || NNTYPE ==  nucleon_t::pn) {
        return &t.pp;
    }
    else if(NNTYPE == nucleon_t::np || NNTYPE ==  nucleon_t::nn){
        return &t.pn;
    }
    else{
        assert(false);
        return nullptr;
    }
}

/*
 * return target part of T type for NNTYPE nucleon-nucleon combination
 * ie.
 * auto t_fe = get_target<nucleon_t::pp>(gm.fe_integrated);
 */
template <nucleon_t NNTYPE, typename T>
const Functional* get_target(const T &t){

    if constexpr (NNTYPE == nucleon_t::pp || NNTYPE ==  nucleon_t::np) {
        return &t.tp;
    }
    else if(NNTYPE == nucleon_t::pn || NNTYPE ==  nucleon_t::nn){
        return &t.tn;
    }
    else {
        assert(false);
        return nullptr;
    }
}

template <typename T>
const Functional* get_target(const T &t, nucleon_t NNTYPE){

    if (NNTYPE == nucleon_t::pp || NNTYPE ==  nucleon_t::np) {
        return &t.tp;
    }
    else if(NNTYPE == nucleon_t::pn || NNTYPE ==  nucleon_t::nn){
        return &t.tn;
    }
    else {
        assert(false);
        return nullptr;
    }
}

/////////////////////////////////////////////////
/// \brief container for phaseshift, 4 components
///
struct phaseshift_4c_type{
    Functional Xpp;
    Functional Xpn;
    Functional Xnp;
    Functional Xnn;
};


/////////////////////////////////////////////////
/// \brief container for phaseshift, 1 components, MOL type
///
struct phaseshift_mol_type{
    Functional _X;
};

/**
 * @tparam F - struct with phase function to be used
 * @param model - glauber model variable which will be calculated
 * @param E - energy for which phase shifts will be calculated
 * calculate phase shift for gm, MOL type
 */
template<typename F, typename GM>
void Calculate_mol(GM &model, double E){
    double Rmax = std::max(model.z_integrated.pp.max(), model.z_integrated.pn.max()) + std::max(model.z_integrated.tp.max(), model.z_integrated.tn.max());
    std::vector<double> x = linspace_vector(0, Rmax, max_b_steps);
    std::vector<double> t;
    t.reserve(max_b_steps);

    for (auto b : x) {
            double s=0;
            if(is_ftype<DiracFunction>(model.z_integrated.tp) || is_ftype<DiracFunction>(model.z_integrated.tn)){
                s+= F::template phase_function<nucleonmol_t::target_p>(model, b, E);
                s+= F::template phase_function<nucleonmol_t::target_n>(model, b, E);
            }
            else{
                s+= F::template phase_function<nucleonmol_t::projectile_p>(model, b, E);
                s+= F::template phase_function<nucleonmol_t::projectile_n>(model, b, E);
            }
            if(is_ftype<DiracFunction>(model.z_integrated.pp) || is_ftype<DiracFunction>(model.z_integrated.pn)){
                s+= F::template phase_function<nucleonmol_t::projectile_p>(model, b, E);
                s+= F::template phase_function<nucleonmol_t::projectile_n>(model, b, E);
            }
            else{
                s+= F::template phase_function<nucleonmol_t::target_p>(model, b, E);
                s+= F::template phase_function<nucleonmol_t::target_n>(model, b, E);
            }

             t.push_back(s);
        }
    model.phase._X = Interpolator(x, t);
}


/**
 * @tparam F - struct with phase function to be used
 * @param model - glauber model variable which will be calculated
 * @param E - energy for which phase shifts will be calculated
 * calculate phase shift for gm, OLA type
 */
template<typename GM>
void Calculate_4c3(GM &model, double E){
    double Rmax = std::max(model.z_integrated.pp.max(), model.z_integrated.pn.max()) + std::max(model.z_integrated.tp.max(), model.z_integrated.tn.max());
    std::vector<double> x = linspace_vector(0, Rmax, max_b_steps);
    std::vector<double> tpp;
    std::vector<double> tnn;
    std::vector<double> tnp;
    std::vector<double> tpn;
    tpp.reserve(max_b_steps);
    tnn.reserve(max_b_steps);
    tnp.reserve(max_b_steps);
    tpn.reserve(max_b_steps);

    #ifdef USE_THREADS
    std::vector<std::thread> threads(4);
    //std::future<void> fut[4]; 
    #endif

    /////////// p-p case
    if (model.projectile.Z() > 1 && model.target.Z() > 1){
        auto fpp = [&]() {
                for (auto b : x)tpp.push_back(model.template phase_function<nucleon_t::pp>(b, E));
                model.phase.Xpp = Interpolator(x, tpp);
                };
        #ifdef USE_THREADS
        threads[0] = std::thread(fpp);
        //fut[0] = std::async(std::launch::async, fpp);
        #else
        fpp();
        #endif
    }
    else if(model.projectile.Z()==0 || model.target.Z()==0){
        model.phase.Xpp = ConstantFunction(0.0);
    }
    else {
        //std::for_each(x.begin(),x.end(),[&](const double &b){
        for (auto b : x)tpp.push_back(model.template phase_function_dirac<nucleon_t::pp>(b, E));
        model.phase.Xpp = Interpolator(x, tpp);
    //});
    }

    //////////////////// p-n case
    if (model.projectile.Z() > 1 && model.target.N() > 1){
        auto fpn = [&]() {
                for (auto b : x)tpn.push_back(model.template phase_function<nucleon_t::pn>(b, E));
                model.phase.Xpn = Interpolator(x, tpn);
                };
        #ifdef USE_THREADS
        threads[1] = std::thread(fpn);
        //fut[1] = std::async(std::launch::async, fpn);
        #else
        fpn();
        #endif
    }
    else if(model.projectile.Z()==0 || model.target.N()==0){
        model.phase.Xpn = ConstantFunction(0.0);
    }
    else {
//        std::for_each(x.begin(),x.end(),[&](const double &b){
        for (auto b : x)tpn.push_back(model.template phase_function_dirac<nucleon_t::pn>(b, E));
        model.phase.Xpn = Interpolator(x, tpn);
  //  });
  }

    ////////////////// n-n case 
    if (model.projectile.N() > 1 && model.target.N() > 1){
        auto fnn  = [&]() {
                for (auto b : x)tnn.push_back(model.template phase_function<nucleon_t::nn>(b, E));
                model.phase.Xnn = Interpolator(x, tnn);
                };
        #ifdef USE_THREADS
        threads[2] = std::thread(fnn);
        //fut[2] = std::async(std::launch::async, fnn);
        #else
        fnn();
        #endif
    }
    else if(model.projectile.N()==0 || model.target.N()==0){
        model.phase.Xnn = ConstantFunction(0.0);
    }
    else {
        //std::for_each(x.begin(),x.end(),[&](const double &b){
        for (auto b : x)tnn.push_back(model.template phase_function_dirac<nucleon_t::nn>(b, E));
        model.phase.Xnn = Interpolator(x, tnn);
    //});
    }

    ////////////////// n-p case
    if (model.projectile.N() > 1 && model.target.Z() > 1){
        auto fnp = [&]() {
                for (auto b : x)tnp.push_back(model.template phase_function<nucleon_t::np>(b, E));
                model.phase.Xnp = Interpolator(x, tnp);
                };
        #ifdef USE_THREADS
        threads[3] = std::thread(fnp);
        //fut[3] = std::async(std::launch::async, fnp);
        #else
        fnp();
        #endif
    }
    else if(model.projectile.N()==0 || model.target.Z()==0){
        model.phase.Xnp = ConstantFunction(0.0);
    }
    else {
        //std::for_each(x.begin(),x.end(),[&](const double &b){
        for (auto b : x)tnp.push_back(model.template phase_function_dirac<nucleon_t::np>(b, E));
        model.phase.Xnp = Interpolator(x, tnp);
    //});
    }

    #ifdef USE_THREADS
        //for(int i=0;i<4;i++){
        //    if(fut[i].valid())fut[i].wait();
        //}
        for (auto& e : threads) {
        if (e.joinable())e.join();
        }
    #endif

} // end of Calculate

}
#endif
