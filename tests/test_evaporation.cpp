#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include <vector>
#include "nurex/ame16.h"
#include "nurex/Utils.h"
#include "nurex/nurex.h"
#include "nurex/evaporation.h"

using namespace nurex;

double energy_effective(int A, int Z, double E){
    prefragment f(A,Z);
    return E - f.pairing;
}

TEST_CASE("W_GS"){
        double n;
        int a = 1;
        auto fx = [&](double x)->double{
            return w_gs(x,a,40);
        };

        auto cdf_fx = [&](double x)->double{
            return cdf_w_gs(x,a,40);
        };

        for(int i=1;i<8;i++){
            a = i;
            n = integrator_adaptive.integrate(fx,0,40*i,1e-5);
            CHECK(n == approx(1,1e-4));
        }

        CHECK(w_gs(0,1,40)==approx(0.05,1e-3));
        CHECK(w_gs(40,1,40)==approx(0.00,1e-3));
        CHECK(w_gs(50,1,40)==approx(0.00,1e-3));
        a = 2;
        CHECK(integrator_adaptive.integrate(fx,0,52,1e-4) == approx(cdf_fx(0),1e-3));
        CHECK(integrator_adaptive.integrate(fx,10,52,1e-4) == approx(cdf_fx(10),1e-3));
        a = 3;
        CHECK(integrator_adaptive.integrate(fx,30,62,1e-4) == approx(cdf_fx(30),1e-3));
        CHECK(w_gs(65,2,40)==0.0);
    }

TEST_CASE("excitation_function"){
    ExcitationFunction exf;
    std::vector<double> e,w;
    e = linspace_vector(0,40,400);
    for(auto en:e){w.push_back(w_gs(en, 1, 38.0));}
    exf.set(e,w,1);
    e = linspace_vector(0,80,100);
    w.clear();
    for(auto en:e){w.push_back(w_gs(en, 2, 38.0));}
    exf.set(e,w,2);

    CHECK(exf.norm[1]==approx(1.0,1e-3));
    CHECK(exf.norm[2]==approx(1.0,1e-3));
    for(auto ee:{0.1, 1.5, 5.5, 11.5, 22.0, 35.0, 38.0}){
        CHECK_MESSAGE(exf.w(ee, 1)==approx(w_gs(ee, 1, 38.),1e-4),"energy:",ee);
        CHECK_MESSAGE(exf.w(ee, 5)==approx(0.0,1e-6),"energy:",ee);
    }

    for(auto ee:{0.1, 11.5, 20.0, 35.0, 50.0, 70.0}){
        CHECK(exf.w(ee, 2)==approx(w_gs(ee, 2, 38.),1e-5));
    }

    // custom excitation integration test
    //Functional fone(ConstantFunction{1.0});
    auto fone = [](double x){return 1.0;};
    auto fhalf = [](double x){return 0.5;};
    custom_excitation.reset();
    
    e = linspace_vector(0,40,400);
    w.clear();
    for(auto en:e){w.push_back(w_gs(en, 1, 40.0));}
    custom_excitation.set(e,w,1);

    e = linspace_vector(0,80,200);
    w.clear();
    for(auto en:e){w.push_back(w_gs(en, 2, 40.0));}
    custom_excitation.set(e,w,2);    

    e = linspace_vector(0,80,200);
    w.clear();
    for(auto en:e){w.push_back(w_gs(en, 3, 40.0));}
    custom_excitation.set(e,w,3);

    e = linspace_vector(0,80,200);
    w.clear();
    for(auto en:e){w.push_back(w_gs(en, 4, 40.0));}
    custom_excitation.set(e,w,4);

    for(auto ee:{0.1, 11.5, 20.0, 35.0, 50.0, 70.0}){                              
        CHECK_MESSAGE(custom_excitation.w(ee, 1)==approx(w_gs(ee, 1, 40.),1e-4),"energy:",ee);
        CHECK_MESSAGE(cdf_wfx_custom(fone, ee, 1)==approx(cdf_wfx_gs(fone, ee, 1, 40.0),1e-3),ee);
        CHECK_MESSAGE(2.0*cdf_wfx_custom(fhalf, ee, 1)==approx(cdf_wfx_gs(fone, ee, 1, 40.0),1e-3),ee);
    }

    // normalisation test
    custom_excitation.reset();
    for(int i=1;i<4;i++){
        std::vector<double> e,w;
        double Ex = 40.0;
        e = linspace_vector(0,80,400);
        for(auto en:e){w.push_back(0.5*w_gs(en, i, Ex));}
        custom_excitation.set(e,w,i);        
        CHECK(custom_excitation.norm[i]==approx(2.0,1e-3));
    }
    

}

TEST_CASE("w_ericson"){
        double n;
        int a = 1;
        auto fx = [&](double x)->double{
            return w_ericson(x,a,40);
        };

        auto cdf_fx = [&](double x)->double{
            return cdf_w_ericson(x,a,40);
        };

        for(int i=1;i<8;i++){
            a = i;
            n = integrator_adaptive.integrate(fx,0,40*i,1e-4);
            CHECK(n == approx(1,1e-4));
        }

        CHECK(w_ericson(0,1,40)==approx(0.025,1e-5));
        CHECK(w_ericson(40,1,40)==approx(0.025,1e-5));
        CHECK(w_ericson(40.01,1,40)==approx(0.00,1e-5));
        a = 2;
        CHECK(integrator_adaptive.integrate(fx,0,80,1e-4) == approx(cdf_fx(0),1e-3));
        CHECK(integrator_adaptive.integrate(fx,10,80,1e-4) == approx(cdf_fx(10),1e-3));
        CHECK(w_ericson(0,2,40)==0.0);
        CHECK(w_ericson(40,2,40)==approx(0.025/2.0,1e-5));
        CHECK(w_ericson(80,2,40)==approx(0.025,1e-5));
        CHECK(w_ericson(80.01,2,40)==approx(0.0,1e-5));
        a = 3;
        CHECK(integrator_adaptive.integrate(fx,30,120,1e-4) == approx(cdf_fx(30),1e-3));

    }

// TEST_CASE("const. temperature approx."){
//     EvaporationParameters c;
//     c.density = level_density_type::ABLA;
//     prefragment f(12,6,c);
//     double r;
//     r = shell_effect_damping(12,asymptotic_level_density_parameter(12, bs_ratio(12,6)),10-f.pairing);
//     CHECK(r == approx(0.3,0.1));

//     // cons. temperature density
//     double a;
//     a = asymptotic_level_density_parameter(12,bs_ratio(12,6));
//     CHECK(a == approx(f.atilda).R(1e-6));
//     r = constant_temperature_parameter(f);
//     CHECK(r == approx(3.1,0.2));

//     a = asymptotic_level_density_parameter(20,bs_ratio(20,6));
//     r = constant_temperature_parameter(prefragment(20,6,c));
//     CHECK(r == approx(2.1,0.1));

//     a = asymptotic_level_density_parameter(50,bs_ratio(50,20));
//     r = constant_temperature_parameter(prefragment(50,20,c));
//     CHECK(r == approx(1.2,0.2));

//     a = asymptotic_level_density_parameter(50,bs_ratio(50,28));
//     r = constant_temperature_parameter(prefragment(50,28,c));
//     CHECK(r == approx(1.19,0.1));

//     a = asymptotic_level_density_parameter(58,bs_ratio(58,28));
//     r = constant_temperature_parameter(prefragment(58,28,c));
//     CHECK(r == approx(1.08,0.1));
//     r = const_temperature_density(f, 10-f.pairing);
//     CHECK(r == approx(0.901342,1e-3));
//     r = const_temperature_density(f, 20-f.pairing);
//     CHECK(r == approx(22.7234,1e-3));
//     r = const_temperature_density(prefragment(132,50,c),50-f.pairing);
//     CHECK(r == approx(1832,10));
// }

// TEST_CASE("fermi gas density"){
//     double r;
//     EvaporationParameters c;
//     c.density = level_density_type::ABLA;
//     prefragment f(12,6,c);
//     ////////// Fermi gas density///////////
//     r = fermi_gas_density(prefragment(12,6,c), 10-f.pairing);
//     CHECK(r == approx(3.4,1e-1));
//     r = fermi_gas_density(f, 29-f.pairing);
//     CHECK(r == approx(394,10));
//     r = fermi_gas_density(f, 30-f.pairing);
//     CHECK(r == approx(472,20));

//     prefragment ca40(40,20,c);
//     r = fermi_gas_density(ca40, 10 - ca40.pairing);
//     CHECK(r == approx(930,100));
    
//     prefragment ca48(48,20,c);
//     r = fermi_gas_density(ca48, 10-ca48.pairing);
//     CHECK(r == approx(1568,150));

//     prefragment sn132(132,50);
//     r = fermi_gas_density(sn132, 12 - sn132.pairing);
//     CHECK(r == approx(200761,150));
// }
/*
TEST_CASE("misc"){
    double r,a;
    prefragment c12(12,6);

    CHECK(J(12,0.0).first == approx(0.808,0.1));
    r = J(18,frdm::get_b2(18,8)).first;
    CHECK(r == approx(1.6,0.1));

    a = asymptotic_level_density_parameter(12,bs_ratio(12,6));
    r = temperature_parameter(energy_corrected(12,6,energy_effective(12,6,10)),a);
    CHECK(r == approx(1.9,0.2));
    r = temperature_parameter(energy_corrected(12,6,energy_effective(12,6,20)),a);
    CHECK(r == approx(3.5,0.5));
    r = temperature_parameter(energy_corrected(12,6,energy_effective(12,6,40)),a);
    CHECK(r == approx(5.2,0.5));

    double e = energy_effective(12,6,10);
    double s2;

    s2 = J(12,0.0).first*constant_temperature_parameter(c12)
        +J(12,0.0).second*constant_temperature_parameter(c12);
    CHECK(s2 == approx(5.0,1.));

    r = J_density_factor(0.0, s2);
    CHECK(r == approx(0.017,0.003));
    r = J_density_factor(1.0, s2);
    CHECK(r == approx(0.043,0.006));

    s2 = J(12,0.0).first*temperature_parameter(energy_corrected(12,6,energy_effective(12,6,40)),a)
        +J(12,0.0).second*temperature_parameter(energy_corrected(12,6,energy_effective(12,6,40)),a);
    CHECK(s2 == approx(8.4,1.));
    r = J_density_factor(1.0, s2);
    CHECK(r == approx(0.02,0.006));

    a = asymptotic_level_density_parameter(132,bs_ratio(132,50));
    s2 = J(132,0.0).first*temperature_parameter(energy_corrected(132,50,energy_effective(132,50,40)),a)
        +J(132,0.0).second*temperature_parameter(energy_corrected(132,50,energy_effective(132,50,40)),a);
    CHECK(s2 == approx(135.,1.));

    // critical energy //
    r = superfluid_phase_critical_energy(prefragment(132, 50), 12);
    CHECK(r == approx(8.,1));

    
    CHECK( S(c12,1,0) == approx(ame16::Sn(12,6),0.001));

    CHECK( S(c12,1,1) == approx(ame16::Sp(12,6),0.001));
    CHECK( C(c12,1,1) == approx(fusion_barrier(12-1,6-1,1,1),0.0001));

    CHECK( S(c12,4,2) == approx(ame16::Sa(12,6),0.001));
    CHECK( C(c12,4,2)== approx(fusion_barrier(12-4,6-2,4,2),0.0001));
}
*/
TEST_CASE("l"){
    auto l = l_orb_distribution(12,4,1,10,asymptotic_level_density_parameter(12));
    CHECK(l.first == approx(4.5,0.5));
    CHECK(l.second == approx(1.7,0.4));

    l = l_orb_distribution(12,1,1,0.27,asymptotic_level_density_parameter(12));
    CHECK(l.first == approx(0.57,0.1));
    CHECK(l.second == approx(0.37,0.1));

    l = l_orb_distribution(12,1,2,0.27,asymptotic_level_density_parameter(12));
    CHECK(l.first == approx(0.95,0.1));
    CHECK(l.second == approx(0.27,0.1));
}
/*
TEST_CASE("level dens."){
    double r;
    nurex::default_evaporation.config = evaporation_config_type::j_single;
    r = level_density(12,6,9,1).first;
    CHECK(r == approx(0.028,0.01));
    r = level_density(12,6,9,0).first;
    CHECK(r == approx(0.011,0.01));
    r = level_density(12,6,29,1).first;
    CHECK(r == approx(10.8,1));
    r = level_density(12,6,29,0).first;
    CHECK(r == approx(4.1,1));
    r = level_density(12,6,40,1).first;
    CHECK(r == approx(54.5,1));
    r = level_density(12,6,40,0).first;
    CHECK(r == approx(20.5,1));

    r = level_density(40,20,9,0).first;
    CHECK(r == approx(0.18,0.05));
    r = level_density(40,20,9,1).first;
    CHECK(r == approx(0.5,0.1));
    r = level_density(40,20,29,1).first;
    CHECK(r == approx(10224,100));
    r = level_density(40,20,29,2).first;
    CHECK(r == approx(16007).R(0.01));
    r = level_density(40,20,29,3).first;
    CHECK(r == approx(20400).R(0.01));
    r = level_density(40,20,29,4).first;
    CHECK(r == approx(23145).R(0.01));
    r = level_density(40,20,29,5).first;
    CHECK(r == approx(24193).R(0.01));
    r = level_density(40,20,29,6).first;
    CHECK(r == approx(23700).R(0.01));
    r = level_density(40,20,29,12).first;
    CHECK(r == approx(7660).R(0.01));

    r = level_density(132,50,12,1).first;
    CHECK(r == approx(66.,3));
    r = level_density(132,50,32,1).first;
    CHECK(r == approx(9e7,1e7));
}
*/
/*
TEST_CASE("width"){
    int A, Z;
    double Ex,j;
    A = 12;
    Z = 6;
    prefragment c12 (A, Z);
    double Sn = ame16::Sn(A, Z);
    double Sa = ame16::Sa(A, Z);
    double Ca = fusion_barrier(A-4, Z-2, 4, 2);
    double Sp = ame16::Sp(A, Z);
    double Cp = fusion_barrier(A-1, Z-1, 1, 1);
    j=0;

    Ex = 11;
    auto g_a = level_density(A-4, Z-2, Ex-Sa-Ca, j);
    double w_a  =   width_e(c12,4,2,Ca,Sa+Ca,Ex, j);
    CHECK(w_a == approx(6.6,0.5));

    Ex = 21;
    auto g_p = level_density(A-1, Z-1, Ex-Sp-Cp, j);
    auto w_p = width_e(c12,1,1,Cp,Sp+Cp,Ex,j);
    CHECK(w_p == approx(5.1,1.0));
    auto g_n = level_density(A-1, Z, Ex-Sn, j);
    auto w_n = width_e(c12,1,1,0.0,Sn,Ex,j);
    CHECK(w_n == approx(4.5,1.0));
    g_a = level_density(A-4, Z-2, Ex-Sa-Ca, j);
    w_a = width_e(c12,4,2,Ca,Sa+Ca,Ex,0);
    CHECK(w_a == approx(13.3,1.0));

    Ex = 29;
    g_p = level_density(A-1, Z-1, Ex-Sp-Cp, j);
    w_p = width_e(c12,1,1,Cp,Sp+Cp,Ex,0);
    CHECK(w_p == approx(6.1,1.0));
    g_n = level_density(A-1, Z, Ex-Sn, j);
    w_n = width_e(c12,1,1,0.0,Sn,Ex);
    CHECK(w_n == approx(8.7,1.0));
    g_a = level_density(A-4, Z-2, Ex-Sa-Ca, j);
    w_a = width_e(c12,4,2,Ca,Sa+Ca,Ex);
    CHECK(w_a == approx(14.1,1.0));

    Ex = 32;
    g_p = level_density(A-1, Z-1, Ex-Sp-Cp, j);
    w_p = width_e(c12,1,1,Cp,Sp+Cp,Ex);
    CHECK(w_p == approx(7.1,1.0));
    g_n = level_density(A-1, Z, Ex-Sn, j);
    w_n = width_e(c12,1,1,0.0,Sn,Ex);
    CHECK(w_n == approx(8.7,1.0));
    g_a = level_density(A-4, Z-2, Ex-Sa-Ca, j);
    w_a = width_e(c12,4,2,Ca,Sa+Ca,Ex);
    CHECK(w_a == approx(16.7,1.0));

    A = 50;
    Z = 20;
    prefragment ca50(A,Z);
    Sn = ame16::Sn(A, Z);
    Sa = ame16::Sa(A, Z);
    Ca = fusion_barrier(A-4, Z-2, 4, 2);
    Sp = ame16::Sp(A, Z);
    Cp = fusion_barrier(A-1, Z-1, 1, 1);

    Ex = 12;
    g_n = level_density(A-1, Z, Ex-Sn, j);
    w_n = width_e(ca50,1,1,0.0,Sn,Ex);
    CHECK(w_n == approx(2.5,0.5));
    g_p = level_density(A-1, Z-1, Ex-Sp-Cp, j);
    w_p = width_e(ca50,1,1,Cp,Sp+Cp,Ex);
    CHECK(w_p == approx(0.0,1.0));
    g_a = level_density(A-4, Z-2, Ex-Sa-Ca, j);
    w_a = width_e(ca50,4,2,Ca,Sa+Ca,Ex);
    CHECK(w_a == approx(0.0,1.0));

    Ex = 22;
    g_n = level_density(A-1, Z, Ex-Sn, j);
    w_n = width_e(ca50,1,1,0.0,Sn,Ex);
    CHECK(w_n == approx(4.8,0.5));
    g_p = level_density(A-1, Z-1, Ex-Sp-Cp, j);
    w_p = width_e(ca50,1,1,Cp,Sp+Cp,Ex);
    CHECK(w_p == approx(0.3,0.2));
    g_a = level_density(A-4, Z-2, Ex-Sa-Ca, j);
    w_a = width_e(ca50,4,2,Ca,Sa+Ca,Ex);
    CHECK(w_a == approx(0.6,0.2));

    Ex = 42;
    g_n = level_density(A-1, Z, Ex-Sn, j);
    w_n = width_e(ca50,1,1,0.0,Sn,Ex);
    CHECK(w_n == approx(9.3,0.5));
    g_p = level_density(A-1, Z-1, Ex-Sp-Cp, j);
    w_p = width_e(ca50,1,1,Cp,Sp+Cp,Ex);
    CHECK(w_p == approx(2.4,0.5));
    g_a = level_density(A-4, Z-2, Ex-Sa-Ca, j);
    w_a = width_e(ca50,4,2,Ca,Sa+Ca,Ex);
    CHECK(w_a == approx(3.3,0.5));

    Ex = 19;
    j = 0   ;
    A = 12;
    Z = 6;
    c12.config.config = evaporation_config_type::j_single;
    c12.config.preset = evaporation_preset_type::abla;
    auto g_CN = level_density(c12, Ex, j).first;

    double r = decay_width(c12,1,0,Ex,0)/g_CN;
    CHECK(r == approx(0.014,0.001));
    r = decay_width(c12,1,0,Ex,0)/g_CN;
    CHECK(r == approx(0.014,0.001));
    r = decay_width(c12,1,1,Ex,0)/g_CN;
    CHECK(r == approx(0.071,0.005));
    r = decay_width(c12,4,2,Ex,0)/g_CN;
    CHECK(r == approx(0.66,0.05));
    r = mean_decay_width(c12,1,0, Ex,0)/g_CN;
    CHECK(r == approx(0.014,0.001));
    r = mean_decay_width(c12,4,2, Ex,0)/g_CN;
    CHECK(r == approx(0.66,0.05));
    Ex = 35;
    g_CN = level_density(c12, Ex, j).first;
    r = decay_width(c12,1,0,Ex,0)/g_CN;
    CHECK(r == approx(0.25,0.05));
    Ex = 55;
    g_CN = level_density(c12, Ex, j).first;
    r = decay_width(c12,1,0,Ex,0)/g_CN;
    CHECK(r == approx(1.8,0.3));

    Ex = 19;
    j = 1;
    g_CN = level_density(c12, Ex,j).first;
    r = decay_width(c12,1,0,Ex,0)/g_CN;
    CHECK(r == approx(0.0056,0.001));
    r = decay_width(c12,1,0,Ex,1)/g_CN;
    CHECK(r == approx(0.000,0.00001));
    r = mean_decay_width(c12,1,0, Ex,j)/g_CN;
    CHECK(r == approx(0.0033,0.0001));

    r = mean_decay_width(c12,2,1, Ex,j);
    r = mean_decay_width(c12,3,1, Ex,j);
    CHECK(r == 0.0);

    // alpha
    r = decay_width(c12,4,2,Ex,0)/g_CN;
    CHECK(r == approx(0.26,0.05));
    r = decay_width(c12,4,2,Ex,1)/g_CN;
    CHECK(r == approx(0.59,0.05));
    r = decay_width(c12,4,2,Ex,2)/g_CN;
    CHECK(r == approx(0.50,0.05));
    r = decay_width(c12,4,2,Ex,3)/g_CN;
    CHECK(r == approx(0.0,0.00001));
    r = decay_width(c12,4,2,Ex,3)/g_CN;
    CHECK(r == approx(0.0,0.00001));
    r = mean_decay_width(c12,4,2,Ex,j)/g_CN;
    CHECK(r == approx(0.13,0.01));
    // proton
    r = mean_decay_width(c12,1,1, Ex,j)/g_CN;
    CHECK(r == approx(0.03,0.01));

    j=2;Ex=29;
    g_CN = level_density(c12, Ex,j).first;
    r = decay_width(c12,1,0,Ex,0)/g_CN;
    CHECK(r == approx(0.03,0.01));
    r = decay_width(c12,1,0,Ex,1)/g_CN;
    CHECK(r == approx(0.07,0.01));
    r = decay_width(c12,1,0,Ex,2)/g_CN;
    CHECK(r == approx(0.077,0.01));
    r = decay_width(c12,1,0,Ex,3)/g_CN;
    CHECK(r == approx(0.046,0.01));
    r = decay_width(c12,1,0,Ex,4)/g_CN;
    CHECK(r == approx(0.00,0.00001));
    r = mean_decay_width(c12,1,0, Ex,j)/g_CN;
    CHECK(r == approx(0.06,0.01));

    // alpha
    r = decay_width(c12,4,2,Ex,0)/g_CN;
    CHECK(r == approx(0.14,0.05));
    r = decay_width(c12,4,2,Ex,1)/g_CN;
    CHECK(r == approx(0.31,0.05));
    r = decay_width(c12,4,2,Ex,2)/g_CN;
    CHECK(r == approx(0.29,0.05));
    r = decay_width(c12,4,2,Ex,3)/g_CN;
    CHECK(r == approx(0.17,0.05));
    r = decay_width(c12,4,2,Ex,4)/g_CN;
    CHECK(r == approx(0.01,0.01));
    r = decay_width(c12,4,2,Ex,5)/g_CN;
    CHECK(r == approx(0.0,0.0001));
    r = mean_decay_width(c12,4,2, Ex,j)/g_CN;
    CHECK(r == approx(0.13,0.02));

    // proton
    r = mean_decay_width(c12,1,1,Ex,j)/g_CN;
    CHECK(r == approx(0.07,0.02));


    //40ca
    A = 40;
    Z = 20;
    prefragment ca40(A,Z);
    j = 2;
    Ex = 20;
    ca40.config.config = evaporation_config_type::j_single;
    ca40.config.preset = evaporation_preset_type::abla;
    g_CN = level_density(ca40, Ex, j).first;
    r = mean_decay_width(ca40,1,1, Ex,j)/g_CN;
    CHECK(r == approx(0.0019,0.0005));
    r = mean_decay_width(ca40,4,2, Ex,j)/g_CN;
    CHECK(r == approx(0.00013,0.00006));
    r = mean_decay_width(ca40,1,0, Ex,j)/g_CN;
    CHECK(r == approx(0.00027,0.0001));

    Ex = 30;
    g_CN = level_density(ca40, Ex, j).first;
    r = mean_decay_width(ca40,1,0, Ex,j)/g_CN;
    CHECK(r == approx(0.03,0.01));
    r = mean_decay_width(ca40,1,1,Ex,j)/g_CN;
    CHECK(r == approx(0.065,0.02));

}
*/
TEST_CASE("P_sanity"){
    int Z = 6;
    for(int A = 9; A<19;A++)
    for(double E = 5; E<150; E = E+2){
        prefragment f1(A,Z);
        f1.config.preset = evaporation_preset_type::nurex;
        auto r1 = evaporation_ratios(f1, E, 0.0);
        CHECK(r1.n.G < 1.0001);
        CHECK(r1.p.G < 1.0001);
        CHECK(r1.g.G < 1.0001);
        CHECK(r1.a.G < 1.0001);
        CHECK(r1.d.G < 1.0001);
        CHECK(r1.he3.G < 1.0001);
        CHECK(r1.imf.G < 1.0001);
    }
}   
// TEST_CASE("width_comparison"){
//     int A = 12;
//     int Z = 6;
//     prefragment f1(A,Z);
//     f1.config.preset = evaporation_preset_type::nurex;
//     prefragment f2(A,Z);
//     f2.config.preset = evaporation_preset_type::abla;

//     for(double E = 25; E<99; E = E+2){
//         auto r1 = evaporation_ratios(f1, E, 0.0);
//         auto r2 = evaporation_ratios(f2, E, 0.0);
//         CHECK(r1.n.G == approx(r2.n.G,0.05));
//     }
// }   

TEST_CASE("gm_evaporation"){
    auto nuc = get_default_nucleus(12,6);
    GlauberModel<OLA> gm(nuc,nuc);

    double r1 = SigmaCC(gm,1000.0);
    gm.SetCCCorrection(cc_correction_t::evaporation);
    double r2 = SigmaCC(gm,1000.0);
    CHECK(r1<r2);

    gm.SetEvaporationParameters({0.0});
    CHECK(gm.evaporation_parameters.Emax == 0.0);
    double r3 = SigmaCC(gm,1000.0);
    CHECK(r3<r2);
    CHECK(r3==approx(r1,0.1));

    gm.SetEvaporationParameters({20.0});
    double r4 = SigmaCC(gm,1000.0);
    CHECK(r4>r1);
};

TEST_CASE("gm_evaporation_custom_ex"){
    auto nuc = get_default_nucleus(12,6);
    GlauberModel<OLA> gm(nuc,nuc);
    EvaporationParameters par1{};
    EvaporationParameters par2{};
    par2.excitation_function = excitation_function_type::CUSTOM;

    gm.SetEvaporationParameters(par2);
    CHECK(gm.evaporation_parameters.excitation_function == excitation_function_type::CUSTOM);
    gm.SetEvaporationParameters(par1);
    CHECK(gm.evaporation_parameters.excitation_function == excitation_function_type::GS);
    gm.SetExcitationFunctionType(excitation_function_type::CUSTOM);
    CHECK(gm.evaporation_parameters.excitation_function == excitation_function_type::CUSTOM);
    gm.SetEvaporationParameters(par1);
    CHECK(gm.evaporation_parameters.excitation_function == excitation_function_type::GS);

    gm.SetEvaporationParameters(par1);
    gm.SetCCCorrection(cc_correction_t::none);
    double r1 = SigmaCC(gm,1000.0);

    gm.SetCCCorrection(cc_correction_t::evaporation);
    double r2 = SigmaCC(gm,1000.0);

    // zero excitation function
    custom_excitation.reset();
    gm.SetExcitationFunctionType(excitation_function_type::CUSTOM);
    gm.SetCCCorrection(cc_correction_t::evaporation);
    double r3 = SigmaCC(gm,1000.0);    
    CHECK(r3 == approx(r1,0.5));

    // reproducing gs by custom exctitation
    for(int i=1;i<4;i++){
        std::vector<double> e,w;
        double Ex = Emax(nuc, par2);
        e = linspace_vector(0,80,400);
        for(auto en:e){w.push_back(w_gs(en, i, Ex));}
        custom_excitation.set(e,w,i); }
    
    double r4 = SigmaCC(gm,1000.0);
    CHECK(r4 > r3);    
    CHECK(r4 == approx(r2).R(5e-3));
    
    // now squeezing GS curve to half range
    custom_excitation.reset();
    for(int i=1;i<4;i++){
        std::vector<double> e,w;
        double Ex = Emax(nuc, par2);
        e = linspace_vector(0,40,400);
        for(auto en:e){w.push_back(w_gs(en*2, i, Ex));}
        custom_excitation.set(e,w,i);        
    }
    double r5 = SigmaCC(gm,1000.0);
    CHECK(r5<r4-3);    
    CHECK(r5>r1);    

};