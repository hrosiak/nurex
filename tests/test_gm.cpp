#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include "nurex/nurex.h"

using namespace nurex;

TEST_CASE("gm_constructor"){
        FermiMotion<NNCrossSectionFit>sigmann_FM(90,90);
        Nucleus carbon(12,6,make_density<DensityHO>(1.548415436,1.6038565),make_density<DensityHO>(1.548415436,1.6038565));
        Nucleus ni(64,28,make_density<DensityFermi>(4.2413494,0.50732378),make_density<DensityFermi>(4.2413494,0.50732378));

        auto gm = make_glauber_model<OLA>(ni, carbon);
        GlauberModel<OLA,FermiMotion<NNCrossSectionFit>> gm2(ni, carbon, sigmann_FM);

        CHECK(gm.projectile.A()==64);
        CHECK(gm.target.A()==12);

        CHECK(gm2.projectile.A()==64);
        CHECK(gm2.target.A()==12);

        GlauberModel<OLA,FermiMotion<NNCrossSectionFit>> gm3(ni, carbon, 0.0);
        CHECK(gm3.range_type() == range_t::ZeroRange);
        CHECK(gm3.get_range().pp == 0.0);
        CHECK(gm3.get_range().pn == 0.0);
        CHECK(gm3.get_range().is_same() == true);
        CHECK(gm3.get_range().is_zero() == true);

        GlauberModel<OLA,FermiMotion<NNCrossSectionFit>> gm4(ni, carbon, 0.39);
        CHECK(gm4.range_type() == range_t::FiniteRange);
        CHECK(gm4.get_range().pp == approx(0.39,1e-6));
        CHECK(gm4.get_range().pn == approx(0.39,1e-6));
        CHECK(gm4.get_range().is_same() == true);
        CHECK(gm4.get_range().is_zero() == false);
        gm4.SetRange(0.0);
        CHECK(gm4.get_range().pp == 0.0);
        CHECK(gm4.get_range().pn == 0.0);
        CHECK(gm4.range_type() == range_t::ZeroRange);
        CHECK(gm4.get_range().is_same() == true);        

        GlauberModel<OLA,FermiMotion<NNCrossSectionFit>> gm5(ni, carbon, range_t::FiniteRange);
        CHECK(gm5.range_type() == range_t::FiniteRange);
        CHECK(gm5.get_range().pp == approx(0.39,1e-6));
        CHECK(gm5.get_range().pn == approx(0.39,1e-6));
        CHECK(gm5.get_range().is_same() == true);        
        CHECK(gm5.range(nucleon_t::pp) == approx(0.39,1e-6));
        CHECK(gm5.range(nucleon_t::pn) == approx(0.39,1e-6));
        gm5.SetRange(0.3, 0.2);
        CHECK(gm5.get_range().pp == approx(0.3,1e-6));
        CHECK(gm5.get_range().pn == approx(0.2,1e-6));
        CHECK(gm5.range(nucleon_t::pp) == approx(0.3,1e-6));
        CHECK(gm5.range(nucleon_t::pn) == approx(0.2,1e-6));
        CHECK(gm5.get_range().is_same() == false);        
        auto r = gm5.GetRange();
        CHECK(r.pp == approx(0.3,0.001));
        CHECK(r.pn == approx(0.2,0.001));
}

TEST_CASE("special case"){
        NNCrossSectionFit sigmann;
        Nucleus p = get_default_nucleus("1H");

        GlauberModelOLA_FM gm1(p, p);
        GlauberModelOLA_FM gm2(p, p, range_t::FiniteRange);
        GlauberModelMOL_FM gm3(p, p, range_t::FiniteRange);

        CHECK(gm1.SigmaR(50.0)==sigmann.pp(50.0));
        CHECK(gm1.SigmaR(150.0)==sigmann.pp(150.0));
        CHECK(gm1.SigmaCC(50.0)==sigmann.pp(50.0));
        CHECK(gm1.SigmaCC(150.0)==sigmann.pp(150.0));
        CHECK(gm2.SigmaR(50.0)==sigmann.pp(50.0));
        CHECK(gm2.SigmaR(150.0)==sigmann.pp(150.0));
        CHECK(gm2.SigmaCC(50.0)==sigmann.pp(50.0));
        CHECK(gm2.SigmaCC(150.0)==sigmann.pp(150.0));
        CHECK(gm3.SigmaR(50.0)==sigmann.pp(50.0));
        CHECK(gm3.SigmaR(150.0)==sigmann.pp(150.0));
        CHECK(gm3.SigmaCC(50.0)==sigmann.pp(50.0));
        CHECK(gm3.SigmaCC(150.0)==sigmann.pp(150.0));
    }
    
TEST_CASE("simple coulomb correction"){
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565),DensityHO(1.548415436,1.6038565));
        Nucleus ni(64,28,DensityFermi(4.24,0.507),DensityFermi(4.24,0.507));
        GlauberModelOLA gm(ni, carbon, range_t::FiniteRange);

        CHECK(coulomb_correction_simple(ni, carbon, 30, 2330.9273/0.9139)==approx(0.9139,0.003));
        CHECK(coulomb_correction_simple(ni, carbon, 100, 1960.77/0.97097)==approx(0.97097,0.001));
        CHECK(coulomb_correction_simple(ni, carbon, 1000, 1820.72/0.9969)==approx(0.9969,0.001));

        double r1 = gm.SigmaR(30);
        CHECK(gm.SigmaR(30)==approx(2392,5.0));
        gm.SetCoulombCorrection(coulomb_correction_t::classic);
        CHECK(gm.SigmaR(30) < r1);
        double r2 = coulomb_correction_simple(ni, carbon, 30, r1)*r1;
        CHECK(gm.SigmaR(30) == approx(r2,0.0001));
    }

TEST_CASE("ola"){
        double r,r3;
        Nucleus carbon(12,6,make_density<DensityHO>(1.548415436,1.6038565),make_density<DensityHO>(1.548415436,1.6038565));
        Nucleus ni(64,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));
        Nucleus p(1,1,DensityDirac(),DensityZero());
        GlauberModelOLA gm(ni, carbon, range_t::FiniteRange);
        GlauberModelOLA gm2(carbon, ni, range_t::FiniteRange);
        GlauberModelOLA gmp1(p, ni, range_t::FiniteRange);
        GlauberModelOLA gmp2(ni, p, range_t::FiniteRange);
        GlauberModel<OLA> gm3(carbon, ni);
        CHECK(gm3.beta()==0.0);
        CHECK(gm3.beta()==0.0);
        gm3.SetRange(0.3901);

        r = gm.SigmaR(30);
        r3 = gm3.SigmaR(30);
        CHECK(r==approx(2392,5.0));
        CHECK(r3==approx(r,1));

        r = gm.SigmaR(100);
        r3 = SigmaR(gm,100);
        CHECK(r==approx(1943,5.0));
        CHECK(r==approx(r3,0.1));

        CHECK(SigmaR(gm3,100)==approx(r,1));

        r = gm.SigmaR(200);
        r3 = SigmaR(gm,200);
        CHECK(r==approx(1805,5.0));
        CHECK(r==approx(r3,0.1));

        // symmetry test
        for(double beta : {0.1,0.2,0.3, 0.4})
            for(double e: {50,100,500,1000}){
            gm.SetRange(beta);
            gm2.SetRange(beta);
            CHECK(gm.SigmaR(e) == approx(gm2.SigmaR(e),1));
            CHECK(gmp1.SigmaR(e) == approx(gmp2.SigmaR(e),1));
        }
    }
    
TEST_CASE("repeat test"){
        double r,r2;
        Nucleus carbon(12,6,make_density<DensityHO>(1.548415436,1.6038565),make_density<DensityHO>(1.548415436,1.6038565));
        Nucleus ni(64,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));

        double e = 30.0;
        double p1, p2;
        for(int i=0;i<10;i++){
            GlauberModelOLA gm(ni, carbon, range_t::FiniteRange);
            GlauberModel<OLA> gm2(carbon, carbon, range_t::FiniteRange);
            r = gm.SigmaR(e);
            r2 = gm2.SigmaR(e);
            CHECK(!std::isnan(r));
            CHECK(!std::isnan(r2));
            CHECK(r==approx(2392,5.0));
        }
    }
    
TEST_CASE("ola range change"){
        Nucleus carbon(12,6,make_density<DensityHO>(1.548415436,1.6038565),make_density<DensityHO>(1.548415436,1.6038565));
        Nucleus ni(64,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));

        GlauberModel<OLA> gm1(carbon, ni);
        GlauberModel<OLA> gm2(carbon, ni, range_t::FiniteRange);
        CHECK(gm1.beta() == 0.0);
        CHECK(gm2.beta() == approx(0.39,1e-9));
        CHECK(gm2.beta() == approx(0.39,1e-9));
        CHECK(SigmaR(gm2, 100.0) > SigmaR(gm1, 100.0));
        CHECK(SigmaR(gm2, 100.0) == approx(1943,5.0));
        gm2.SetRange(0.0);
        CHECK(gm2.beta() == 0.0);
        CHECK(gm2.beta() == 0.0);
        CHECK(SigmaR(gm2, 100.0) == approx(SigmaR(gm1, 100.0),0.1));
        gm1.SetRange(0.39);
        CHECK(gm1.beta() == approx(0.39,1e-9));
        CHECK(gm1.beta() == approx(0.39,1e-9));
        CHECK(SigmaR(gm1, 100.0) == approx(1943,5.0));
    }
    
TEST_CASE("ola_fm"){
        double r,r2,r3,r4;
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565),DensityHO(1.548415436,1.6038565));
        Nucleus ni(64,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));
        Nucleus p(1,1,DensityDirac(),DensityZero());
        //64ni on 12C case
        GlauberModelOLA_FM gm(ni, carbon, range_t::FiniteRange);
        GlauberModelOLA_FM gm2(carbon, ni, range_t::FiniteRange);
        GlauberModelOLA_FM gmp1(p, ni, range_t::FiniteRange);
        GlauberModelOLA_FM gmp2(ni, p, range_t::FiniteRange);

        r = gm.SigmaR(30);
        r2 = gm2.SigmaR(30);
        CHECK(r==approx(2745,5.0));
        CHECK(r==approx(r2,1));

        r = gm.SigmaR(100);
        r2 = gm2.SigmaR(100);
        CHECK(r==approx(2106,5.0));
        CHECK(r==approx(r2,1));
        r3 = SigmaR(gm,100);
        r4 = SigmaR(gm2,100);
        CHECK(r==approx(r3,0.1));
        CHECK(r2==approx(r4,0.1));

        r = gm.SigmaR(200);
        r2 = gm2.SigmaR(200);
        CHECK(r==approx(1835,5.0));
        CHECK(r==approx(r2,1));
        r3 = SigmaR(gm,200);
        r4 = SigmaR(gm2,200);
        CHECK(r==approx(r3,0.1));
        CHECK(r2==approx(r4,0.1));

        //12C on 12C case
        GlauberModelOLA_FM gmc(carbon, carbon,range_t::FiniteRange);
        CHECK(gmc.SigmaR(30) == approx(1401,5.0));
        CHECK(gmc.SigmaR(100) == approx(1010,5.0));
        CHECK(gmc.SigmaR(200) == approx(837,5.0));

        r3 = SigmaR(gmc,200);
        CHECK(SigmaR(gmc,200)==approx(gmc.SigmaR(200),0.1));

        // symmetry test
        for(double beta : {0.0,0.1,0.2,0.3, 0.4})
            for(double e: {50,100,500,1000}){
            gm.SetRange(beta);
            gm2.SetRange(beta);
            CHECK(gm.SigmaR(e) == approx(gm2.SigmaR(e),1)); //NN
            CHECK(gmp1.SigmaR(e) == approx(gmp2.SigmaR(e),1)); //pN
        }


        // charge changing 12c on 12c
        CHECK(SigmaCC(gmc,30)==approx(1265.9,6));
        CHECK(SigmaCC(gmc,100)==approx(867.0,4));
        CHECK(SigmaCC(gmc,200)==approx(689.7,3));
    }

TEST_CASE("ola_fmd"){
        const double rel_prec = 0.0075;
        double r,r2,r3,r4;
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565),DensityHO(1.548415436,1.6038565));
        Nucleus ni(64,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));
        Nucleus p(1,1,DensityDirac(),DensityZero());
        
        //64ni on 12C case
        GlauberModel<OLA_FMD, NNCrossSection_FermiMotionD> gm(ni, carbon, range_t::FiniteRange);
        GlauberModel<OLA_FMD, NNCrossSection_FermiMotionD> gm2(carbon, ni, range_t::FiniteRange);
        GlauberModelOLA gm3(ni, carbon, range_t::FiniteRange);
        // proton case
        GlauberModel<OLA_FMD, NNCrossSection_FermiMotionD> gmp1(p, ni, range_t::FiniteRange);
        GlauberModel<OLA_FMD, NNCrossSection_FermiMotionD> gmp2(ni, p, range_t::FiniteRange);

        gm.SetFECoefficient(0.05);
        gm2.SetFECoefficient(0.05);

        r = gm.SigmaR(30);
        r2 = gm2.SigmaR(30);
        CHECK(r==approx(2392).R(rel_prec));
        CHECK(r==approx(r2).R(rel_prec));

        r = gm.SigmaR(100);
        r2 = gm2.SigmaR(100);
        CHECK(r==approx(1943).R(rel_prec));
        CHECK(r==approx(r2).R(rel_prec));

        for(double e : {50,150,300,500}){
            double r1 = gm.SigmaR(e);
            double r2 = gm2.SigmaR(e);
            double r3 = gm3.SigmaR(e);
            CHECK(r1 == approx(r3).R(rel_prec));
            CHECK(r1 == approx(r2).R(rel_prec));
        }

        gm.SetFECoefficient(0.65);
        gm2.SetFECoefficient(0.65);
        gm.SetRange(0.39);
        gm2.SetRange(0.39);

        // check sigma_nn from density symmetry         
        for(double e : {50,150})
            for(double r1:{0.,0.5,1.0,2.0}){
                    CHECK(gm.template sigma_nn_density<nucleon_t::pp>(r1,r1,e) == approx(gm2.template sigma_nn_density<nucleon_t::pp>(r1,r1,e)).R(0.001) );
                    CHECK(gm.template sigma_nn_density<nucleon_t::nn>(r1,r1,e) == approx(gm2.template sigma_nn_density<nucleon_t::nn>(r1,r1,e)).R(0.001) );
                    CHECK(gm.template sigma_nn_density<nucleon_t::pn>(r1,r1,e) == approx(gm2.template sigma_nn_density<nucleon_t::np>(r1,r1,e)).R(0.001) );
                    CHECK(gm.template sigma_nn_density<nucleon_t::np>(r1,r1,e) == approx(gm2.template sigma_nn_density<nucleon_t::pn>(r1,r1,e)).R(0.001) );  
            }
        
        for(double e : {50,150,300,500}){
            double r1 = gm.SigmaR(e);
            double r2 = gm2.SigmaR(e);
            CHECK(r1 == approx(r2,1).R(rel_prec));
            double r3 = gmp1.SigmaR(e);
            double r4 = gmp2.SigmaR(e);
            CHECK(r3 == approx(r4,1).R(rel_prec));
        }
    }
    
TEST_CASE("mol"){
        Nucleus carbon(12,6,make_density<DensityHO>(1.548415436,1.6038565),make_density<DensityHO>(1.548415436,1.6038565));
        Nucleus ni64(64,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));
        Nucleus p(1,1,DensityDirac(),DensityZero());
        
        GlauberModel<MOL> gm(carbon, ni64, range_t::FiniteRange);
        GlauberModel<MOL> gm2(ni64, carbon, range_t::FiniteRange);
        GlauberModel<MOL> gmp1(p, ni64, range_t::FiniteRange);
        GlauberModel<MOL> gmp2(ni64, p, range_t::FiniteRange);

        CHECK(gm.SigmaR(30) == approx(2265.5,4));
        CHECK(gm.SigmaR(100) == approx(1876,3));
        CHECK(gm.SigmaR(200) == approx(1752,2));

        CHECK(gm2.SigmaR(30) == approx(2265.5,4));
        CHECK(gm2.SigmaR(100) == approx(1876,3));
        CHECK(gm2.SigmaR(200) == approx(1752,2));

        CHECK(gm.SigmaCC(30) == approx(-1.0,1e-6)); // MOL does not calculate CC at the moment
        
        const double rel_prec = 0.005;
        for(double e : {30, 50,150,300,500}){
            CHECK(gm.SigmaR(30) == approx(gm2.SigmaR(30)).R(rel_prec));
            CHECK(gmp1.SigmaR(30) == approx(gmp2.SigmaR(30)).R(rel_prec));
        }
        gm.SetRange(0.0);
        gm2.SetRange(0.0);
        gmp1.SetRange(0.0);
        gmp2.SetRange(0.0);
        for(double e : {30, 50,150,300,500}){
            CHECK(gm.SigmaR(30) == approx(gm2.SigmaR(30)).R(rel_prec));
            CHECK(gmp1.SigmaR(30) == approx(gmp2.SigmaR(30)).R(rel_prec));
        }    
    }
    
TEST_CASE("mol4c"){
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565), DensityHO(1.548415436,1.6038565));
        Nucleus ni(64,28,DensityFermi(4.2413494,0.50732378),DensityFermi(4.2413494,0.50732378));
        
        GlauberModelMOL4C gm(ni, carbon, range_t::FiniteRange);
        GlauberModelMOL4C gms(carbon, ni, range_t::FiniteRange);
        GlauberModel<MOL4C> gm2(ni, carbon, range_t::FiniteRange);
        GlauberModel<MOL4C> gm3(carbon,ni, range_t::FiniteRange);

        double r,r2;

        //64Ni on 12C cases
        r = gm.SigmaR(30);
        r2 = gms.SigmaR(30);
        CHECK(r==approx(2304).epsilon(5));
        CHECK(r==approx(r2).epsilon(0.1));
        CHECK(gm2.SigmaR(30.0)==approx(r2).epsilon(1));

        r= gm.SigmaR(100);
        r2 = gms.SigmaR(100);
        CHECK(r==approx(1904).epsilon(5));
        CHECK(r==approx(r2).epsilon(0.1));
        CHECK(gm2.SigmaR(100.0)==approx(r2).epsilon(1));

        r= gm.SigmaR(200);
        r2 = gms.SigmaR(200);
        CHECK(r==approx(1778).epsilon(5));
        CHECK(r==approx(r2).epsilon(0.1));
        CHECK(gm2.SigmaR(200.0)==approx(r2).epsilon(1));

        CHECK(gm.SigmaR(30) == approx(gm3.SigmaR(30),1));
        CHECK(gm.SigmaR(100) == approx(gm3.SigmaR(100),1));
        CHECK(gm.SigmaR(200) == approx(gm3.SigmaR(200),1));
    }
    
TEST_CASE("MOL4C_FM"){
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565), DensityHO(1.548415436,1.6038565));
        Nucleus ni(64,28, DensityFermi(4.2413494,0.50732378), DensityFermi(4.2413494,0.50732378));
        Nucleus p(1,1,DensityDirac(),DensityZero());
        
        GlauberModelMOL4C_FM gm(ni, carbon, range_t::FiniteRange);
        GlauberModelMOL4C_FM gms(carbon, ni, range_t::FiniteRange);
        GlauberModelMOL4C_FM gmc(carbon, carbon, range_t::FiniteRange);
        GlauberModelMOL4C_FM gmp1(ni, p, range_t::FiniteRange);
        GlauberModelMOL4C_FM gmp2(p, ni, range_t::FiniteRange);


        double r,r3,r4;

        //12C on 12C case
        CHECK(gmc.SigmaR(30)==approx(1307).epsilon(1));
        CHECK(gmc.SigmaR(100)==approx(974).epsilon(1));
        CHECK(gmc.SigmaR(200)==approx(817).epsilon(1));

        //64Ni on 12C cases
        CHECK(gm.SigmaR(30)==approx(2600).epsilon(10));
        CHECK(gm.SigmaR(100)==approx(2050).epsilon(10));
        CHECK(gm.SigmaR(200)==approx(1806).epsilon(10));

        //symmetry test
        const double rel_prec = 0.005;
        for(double e : {30, 50,150,300,500}){
            CHECK(gm.SigmaR(30) == approx(gms.SigmaR(30)).R(rel_prec));
            CHECK(gmp1.SigmaR(30) == approx(gmp2.SigmaR(30)).R(rel_prec));
        }
        gm.SetRange(0.0);
        gms.SetRange(0.0);
        gmp1.SetRange(0.0);
        gmp2.SetRange(0.0);
        for(double e : {30, 50,150,300,500}){
            CHECK(gm.SigmaR(30) == approx(gms.SigmaR(30)).R(rel_prec));
            CHECK(gmp1.SigmaR(30) == approx(gmp2.SigmaR(30)).R(rel_prec));
        }    
    }
    
TEST_CASE("pN"){
    Nucleus h = get_default_nucleus(1,1);
    Nucleus c = get_default_nucleus(12,6);
    Nucleus sn = get_default_nucleus(112,50);

    GlauberModel<OLA> gm1(h,c, range_t::FiniteRange);
    GlauberModel<OLA> gm1r(c,h, range_t::FiniteRange);

    CHECK(gm1.SigmaR(30)==approx(421).epsilon(10));
    CHECK(gm1.SigmaR(50)==approx(353).epsilon(10));
    CHECK(gm1.SigmaR(70)==approx(309).epsilon(10));
    CHECK(gm1.SigmaR(100)==approx(267).epsilon(10));
    CHECK(gm1.SigmaR(200)==approx(217).epsilon(10));
    CHECK(gm1.SigmaR(300)==approx(204).epsilon(10));
    CHECK(gm1.SigmaR(500)==approx(212).epsilon(10));
    CHECK(gm1.SigmaR(650)==approx(237).epsilon(10));
    CHECK(gm1.SigmaR(750)==approx(244).epsilon(10));
    CHECK(gm1.SigmaR(1000)==approx(246.5).epsilon(10));
    CHECK(gm1.SigmaR(2000)==approx(247).epsilon(10));
    CHECK(gm1.SigmaR(1000)==approx(gm1r.SigmaR(1000)).epsilon(0.1));

    GlauberModel<OLA> gm2a(h,sn, range_t::FiniteRange);
    GlauberModel<OLA> gm2b(sn,h, range_t::FiniteRange);
    gm2a.SetCoulombCorrection(coulomb_correction_t::none);
    gm2b.SetCoulombCorrection(coulomb_correction_t::none);
    for(double e:{50,100,200,300,500,750,1000}){
        double r1, r2;
        r1 = gm2a.SigmaR(e);
        r2 = gm2b.SigmaR(e);
        CHECK(r1==approx(r2).epsilon(0.1));
    }

    gm2a.SetCoulombCorrection(coulomb_correction_t::classic);
    gm2b.SetCoulombCorrection(coulomb_correction_t::classic);
    for(double e:{50,100,200,300,500,750,1000}){
        CHECK(gm2a.SigmaR(1000)==approx(gm2b.SigmaR(1000)).epsilon(0.1));
    }
    }
    
TEST_CASE("Compound"){
        Nucleus c = get_default_nucleus(12,6);
        Nucleus h = get_default_nucleus(1,1);
        Compound ch2 = Compound({{c,1},{h,2}});

        GlauberModelOLA gmc(c,c,range_t::FiniteRange);
        GlauberModelOLA gmh(c,h,range_t::FiniteRange);
        GMCompound<GlauberModelOLA> gm(c,ch2, range_t::FiniteRange);

        double rc = SigmaR(gmc,100);
        double rh = SigmaR(gmh,100);
        CHECK(rc>rh);

        CHECK(gm.num()==2);
        CHECK(gm.projectile == c);
        CHECK(gm.target.get_nucleus(0)==c);
        CHECK(gm.target.get_nucleus(1)==h);
        CHECK(gm.target.get_nucleus(0).A()==12);
        CHECK(gm.target.get_nucleus(0).Z()==6);
        CHECK(gm.target.get_nucleus(1).A()==1);
        CHECK(gm.target.molar_fraction(0)==1);
        CHECK(gm.target.molar_fraction(1)==2);

        CHECK( (gm.gm[0].projectile)==c);
        CHECK( (gm.gm[1].projectile)==c);
        CHECK( gm.gm[1].projectile.A() == 12);
        CHECK( gm.gm[0].target.A()==12);
        CHECK( gm.gm[0].target.Z()==6);
        //CHECK( *(gm.gm[0].target)==c);
        //CHECK( *(gm.gm[1].target)==h);

        double r1 = SigmaR(gm,100);
        double r2 = rc + 2*rh;
        CHECK(r1 == approx(r2/3.0).R(1e-4));
}

TEST_CASE("sigma_cc + simple correction"){
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565),DensityHO(1.548415436,1.6038565));
        GlauberModelOLA_FM gmc(carbon, carbon, range_t::FiniteRange);
        CHECK(SigmaCC(gmc,30)==approx(1265.9,6));
        CHECK(SigmaCC(gmc,100)==approx(867.0,4));
        CHECK(SigmaCC(gmc,200)==approx(689.7,3));

        gmc.SetCCCorrection(cc_correction_t::PRC82);
        CHECK(SigmaCC(gmc,30)==approx(1265.9*sigma_cc_scaling_factor(30),6));
        CHECK(SigmaCC(gmc,100)==approx(867.0*sigma_cc_scaling_factor(100),4));
        CHECK(SigmaCC(gmc,200)==approx(689.7*sigma_cc_scaling_factor(200),3));
}

TEST_CASE("fmd"){
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565),DensityHO(1.548415436,1.6038565));
        GlauberModelOLA_FM gm1(carbon, carbon, range_t::FiniteRange);
        GlauberModelOLA gm2(carbon, carbon, range_t::FiniteRange);
        GlauberModel<OLA_FMD, FermiMotionD<NNCrossSectionFit>> gm(carbon, carbon, range_t::FiniteRange);
        gm.fe_coefficient = 0.5;
 
        for(double e: {50, 200, 400}){
            CHECK(SigmaR(gm1,e)>SigmaR(gm,e));
            CHECK(SigmaR(gm2,e)<SigmaR(gm,e));
            }

        // should be small effect at high energies
        double e = 1500;
        CHECK(SigmaR(gm1,e) == approx(SigmaR(gm,e),5));
        CHECK(SigmaR(gm2,e) == approx(SigmaR(gm,e),5));
    }
    
TEST_CASE("sigma_xn ola_fm"){
        double r,r2,r3,r4;
        Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565),DensityHO(1.548415436,1.6038565));
        GlauberModelOLA_FM gmc(carbon, carbon,range_t::FiniteRange);
        // charge changing 12c on 12c
        CHECK(SigmaCC(gmc,30)==approx(1265.9,6));
        CHECK(SigmaCC(gmc,100)==approx(867.0,4));
        CHECK(SigmaCC(gmc,200)==approx(689.7,3));

        for(double e:{300,100,200}){
            double r1 = gmc.SigmaR(e);
            double r2 = gmc.SigmaCC(e);
            double r3 = gmc.SigmaXN(e);
            CHECK(r3 == approx(r1-r2,0.1));
            CHECK(SigmaXN(gmc,e)==approx(gmc.SigmaXN(e),0.1));
        }

       auto c12 = get_default_nucleus(12,6);
       auto ca48 = get_default_nucleus(48,20);

       GlauberModelOLA_FM gm(c12, c12);
       r = SigmaIN(gm,300,1);
       CHECK(r == approx(103,5));
       r = SigmaIN(gm,300,2);
       CHECK(r == approx(26,5));

       GlauberModelOLA_FM gm2(ca48, c12);
       r = SigmaXN(gm2,300);
       CHECK(r == approx(273,3));
       r = SigmaIN(gm2,300,2);
       CHECK(r == approx(56,5));
       
       auto rs = SigmaINs(gm,300);
       auto rs_total = SigmaXN(gm,300);
       double sum = 0;
       for(int i=0; i< rs.size();i++){
           double nn = SigmaIN(gm,300,i+1);
           sum += nn;
           CHECK( rs[i] == approx(nn,0.001));
       }
       CHECK(sum == approx(rs_total,0.1));
    }

TEST_CASE("n removal scaling"){
       double r,r2,r3,r4;
       Nucleus carbon(12,6,DensityHO(1.548415436,1.6038565),DensityHO(1.548415436,1.6038565));
       GlauberModelOLA_FM gmc(carbon, carbon,range_t::FiniteRange);
        
       auto c12 = get_default_nucleus(12,6);
       auto ca48 = get_default_nucleus(48,20);
       GlauberModelOLA_FM gm(c12, c12);
       GlauberModelOLA_FM gm2(ca48, c12);
       EvaporationParameters par0;       
       EvaporationParameters par;       
       par.n_removal_scaling = 0.75;
            
       gm.SetEvaporationParameters(par0);
       gm2.SetEvaporationParameters(par0);
       auto rs = SigmaINs(gm,300);
       auto rs2 = SigmaINs(gm2,300);
       double xn0 = SigmaXN(gm, 100);
       double xn2 = SigmaXN(gm2, 100);
       r = SigmaIN(gm,300,1);
       CHECK(gm.evaporation_parameters.n_removal_scaling == 1.0);
        
       gm.SetEvaporationParameters(par);
       gm2.SetEvaporationParameters(par);
       CHECK(gm.evaporation_parameters.n_removal_scaling == 0.75);
       auto rs3 = SigmaINs(gm,300);
       auto rs4 = SigmaINs(gm2,300); 
       double xn1 = SigmaXN(gm, 100); 
       double xn4 = SigmaXN(gm2, 100);            
       r2 = SigmaIN(gm,300,1);

       for(int i=0; i< rs.size();i++){       
           CHECK( 0.75*rs[i] == approx(rs3[i],0.001));
           CHECK( 0.75*rs2[i] == approx(rs4[i],0.001));
       }

        CHECK(0.75*r == approx(r2,0.001));

       CHECK(0.75*xn0 == approx(xn1,0.001));
       CHECK(0.75*xn2 == approx(xn4,0.001));
        
       gm.SetEvaporationParameters(par0);
       gm2.SetEvaporationParameters(par0);
       double cc0 = SigmaCC(gm,100);
       double cc2 = SigmaCC(gm2,100);
       gm.SetCCCorrection(cc_correction_t::evaporation);
       gm2.SetCCCorrection(cc_correction_t::evaporation);
       double cc1 = SigmaCC(gm,100);
       double cc3 = SigmaCC(gm2,100);
       gm.SetEvaporationParameters(par);
       gm2.SetEvaporationParameters(par);
       
       double cc4 = SigmaCC(gm,100);
       double cc5 = SigmaCC(gm2,100);

       CHECK(cc0 < cc1);
       CHECK(cc2 < cc3);
       CHECK(cc4 < cc1);
       CHECK(cc5 < cc3);
       CHECK(cc0 < cc4);
       CHECK(cc2 < cc5);       
    }


