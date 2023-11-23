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
#include "nurex/data_frdm.h"

using namespace nurex;

TEST_CASE("basic"){
        auto fx = [&](double x)->double{
            return x*x*x - x -2;
        };
        double r = bisection(fx,-2,4);
        CHECK(r==approx(1.521,0.001));
        CHECK(p_from_T(1000.0,112.0)==approx(189300,1000));
        CHECK(T_from_p(189300,112.0)==approx(1000,10));

        for(double e:{50,100,200,300,500,1000,2000}){
            double g, b, p, p2, t, e2;
            g = gamma_from_T(e);
            b =  beta_from_gamma(g);
            p = p_from_beta(b);
            p2 = p_from_T(e);
            e2 = T_from_p(p);
            CHECK(e == approx(e2).R(1e-6));
            CHECK(p == approx(p2).R(1e-6));
            p = p_from_beta(b,100);
            p2 = p_from_T(e,100);
            e2 = T_from_p(p,100.0);
            CHECK(p == approx(p2).R(1e-6));
            CHECK(e == approx(e2).R(1e-6));
        }

    }

TEST_CASE("ame"){
        CHECK(ame16::get_mass(12,6)==approx(12.0).R(1e-5));
        CHECK(ame16::get_mass(11,6)==approx(11.011432591,1e-6));
        CHECK(ame16::get_mass(153,63)==approx(152.922,1e-3));
        CHECK(ame16::get_mass(153,12)==0.0);

        CHECK(ame16::get_nuclear_mass(12,6)==approx(11.99671).epsilon(1e-4));
        CHECK(ame16::get_nuclear_mass(11,6)==approx(11.00814).epsilon(1e-4));
        CHECK(ame16::get_nuclear_mass(1,1)==approx(1.00728).R(1e-5));
        CHECK(ame16::get_nuclear_mass(10,1)==0.0);

        CHECK(ame16::Sn(12,6)==approx(18.721,1e-2));
        CHECK(ame16::Sp(12,6)==approx(15.957,1e-2));
        CHECK(ame16::S(12,6,0,1)== approx(18.721,1e-2));
        CHECK(ame16::S(12,6,1,0)==approx(15.957,1e-2));

        CHECK(ame16::Sn(12,6,2)==approx(31.841,1e-2));
        CHECK(ame16::Sp(12,6,2)==approx(27.185,1e-2));
        CHECK(ame16::S(12,6,0,2)==approx(31.841,1e-2));
        CHECK(ame16::S(12,6,2,0)==approx(27.185,1e-2));

        CHECK(ame16::Sn(11,6,2)== approx(34.404,1e-2));
        CHECK(ame16::Sp(11,6,2)== approx(15.277,1e-2));

        CHECK(ame16::Sn(15,6)==approx(1.218,1e-2));
        CHECK(ame16::Sp(15,6)==approx(21.080,1e-2));
        CHECK(ame16::Sn(15,6,2)==approx(9.395,1e-2));
        CHECK(ame16::Sp(15,6,2)==approx(38.364,1e-2));

        CHECK(ame16::Sn(9,6)==approx(14.225,1e-2));
        CHECK(ame16::Sp(9,6)==approx(1.300,1e-2));
        CHECK(ame16::Sn(9,6,2)==approx(0.0,1e-5));
        CHECK(ame16::Sp(9,6,2)==approx(1.436,1e-2));
        CHECK(ame16::S(9,6,0,1)==approx(14.225,1e-2));
        CHECK(ame16::S(9,6,1,0)==approx(1.300,1e-2));
        CHECK(ame16::S(9,6,0,2)==approx(0.0,1e-5));
        CHECK(ame16::S(9,6,2,0)==approx(1.436,1e-2));

        CHECK(ame16::Sn(8,6)== 0.0);
        CHECK(ame16::Sp(8,6)==approx(-0.099,1e-1));
        CHECK(ame16::Sn(8,6,2)== 0.0);
        CHECK(ame16::Sp(8,6,2)==approx(-2.111,1e-2));

        CHECK(ame16::Sn(7,6)== 0.0);
        CHECK(ame16::Sp(7,6)== 0.0);
        CHECK(ame16::Sn(7,6,2)== 0.0);
        CHECK(ame16::Sp(7,6,2)== 0.0);
        CHECK(ame16::S(7,6,2,0)== 0.0);

        CHECK(ame16::Sa(11,6)== approx(7.544,1e-2));
        CHECK(ame16::Sa(10,6)== approx(5.101,1e-2));
        CHECK(ame16::S(10,6,2,2)== approx(5.101,1e-2));
        CHECK(ame16::S(11,6,2,2)== approx(7.544,1e-2));

        CHECK(ame16::BA(1,1)== 0.0);
        CHECK(ame16::BA(2,1)== approx(1.112,1e-2));
        CHECK(ame16::BA(10,6)== approx(6.032,1e-2));
        CHECK(ame16::BA(9,6)== approx(4.337,1e-2));
        CHECK(ame16::BA(8,6)== approx(3.101,1e-2));
        CHECK(ame16::BA(7,6)== 0.0);
    }

TEST_CASE("fermi eneryg"){
        double val = fermi_energy(0.085,proton_mass);
        CHECK(val==approx(38,1));
        CHECK(fermi_momentum(0.085)==approx(sqrt(val*2*proton_mass)).R(1e-4));
    }

TEST_CASE("parametrization"){
        auto c = get_default_nucleus(12,6);
        auto al = get_default_nucleus(27,13);
        double r1, r2;
        r1 = SigmaR_Kox(12,6,30,12,6);
        //CHECK(r1 == approx(1295,20),"%lf vs 1295 +-20 ", r1, r2); // TODO: check this numbers
        CHECK(SigmaR_Kox(12,6,30,12,6) == SigmaR_Kox(c,30,c));
        CHECK(SigmaR_Kox(12,6,83,12,6) == approx(957,20));
        CHECK(SigmaR_Kox(12,6,200,12,6) == approx(890,20));
        CHECK(SigmaR_Kox(12,6,870,12,6) == approx(846,20));
        CHECK(SigmaR_Kox(12,6,870,12,6) == SigmaR_Kox(c,870,c));

        //CHECK(SigmaR_Kox(20,10,20,27,13) == approx(2098,30)); // TODO: check this numbers
        CHECK(SigmaR_Kox(20,10,100,27,13) == approx(1684,30));
        CHECK(SigmaR_Kox(20,10,300,27,13) == approx(1612,30));

        CHECK(SigmaR_Kox(12,6,870,27,13) == SigmaR_Kox(c,870,al));
    }

/*
TEST_CASE("coul. corr."){
    Nucleus sn = get_default_nucleus(112,50);
    Nucleus p = get_default_nucleus(1,1);

    double v1, v2;
    double b = beta_from_T(30);
    v1 = 0.5*closest_distance(112,50,50,b);
    v2 = sommerfeld(50,50,Ecm_from_T(30,112,112));
    CHECK(v1 == approx(v2).R(1e-2));
    v1 = 0.5*closest_distance(112,50,50,b);
    v2 = sommerfeld(50,50,Ecm_from_T(30,112,112));
    CHECK(v1 == approx(v2).R(1e-2));
    b = beta_from_T(60);
    v1 = 0.5*closest_distance(112,50,50,b);
    v2 = sommerfeld(50,50,Ecm_from_T(60,112,112));
    CHECK(v1 == approx(v2).R(1e-2));
    v1 = 0.5*closest_distance(112,50,50,b);
    v2 = sommerfeld(50,50,Ecm_from_T(60,112,112));
    CHECK(v1 == approx(v2).R(1e-2));

    b = beta_from_T(30);
    v1 = 0.5*closest_distance(112,50,1,b);
    v2 = sommerfeld(50,1,Ecm_from_T(30,112,1));
    CHECK(v1 == approx(v2).R(1e-2));
    v1 = 0.5*closest_distance(50,112,1,b);
    v2 = sommerfeld(1,50,Ecm_from_T(30,1,112));
    CHECK(v1 == approx(v2).R(1e-2));

    for(double e:{50,100,200,500,750,1000})
        for(double cs:{200,500,750,1000,1100,1500,2000,2500}){
            double r1, r2;
            r1 = Ecm_from_T_relativistic(e, p.A(), sn.A());
            r2 = Ecm_from_T_relativistic(e, sn.A(), p.A());
            CHECK(r1 == approx(r2,0.001));
            r1 = coulomb_correction_relativistic(sn,p,e, cs);
            r2 = coulomb_correction_relativistic(p,sn, e, cs);
            CHECK(r1 == approx(r2,0.001));
        }

    for(double e:{30,50,70})
        for(double cs:{200,500,750,1000,1100,1500,2000,2500}){
            double r1, r2, r3;
            r1 = Ecm_from_T_relativistic(e, p.A(), sn.A());
            r2 = Ecm_from_T_relativistic(e, sn.A(), p.A());
            r3 = Ecm_from_T(e, sn.A(), p.A());
            CHECK(r1 == approx(r3,1));
            CHECK(r2 == approx(r3,1));
    }
}
*/
TEST_CASE("frdm"){
    using frdm::get_b2, frdm::get_b4, frdm::get_data;

    CHECK(get_b2(16,8)==approx(-0.01,1e-3));
    CHECK(get_b2(53,20)==approx(0.064,1e-3));
    CHECK(get_b2(500,92)==0.0);
    CHECK(get_b2(12,6)==0.0);
    CHECK(get_b4(16,8)==approx(-0.122,1e-4));
    CHECK(get_b4(1,1)==approx(0.0,1e-5));
    CHECK(get_b4(145,47)==approx(0.054,1e-4));
    
    CHECK(get_data(16,8).b2 == approx(-0.01,1e-3));
    CHECK(get_data(53,20).b2==approx(0.064,1e-3));
    CHECK(get_data(500,92).b2==0.0);
    CHECK(get_data(12,6).b2==0.0);
    CHECK(get_data(16,8).b4==approx(-0.122,1e-4));
    CHECK(get_data(1,1).b4==approx(0.0,1e-5));
    CHECK(get_data(145,47).b4==approx(0.054,1e-4));
}

