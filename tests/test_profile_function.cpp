#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include "nurex/nurex.h"

using namespace nurex;


TEST_CASE("OLA"){
        const double rel_prec = 0.001;
        Nucleus projectile = get_default_nucleus(12,6);
        Nucleus target = get_default_nucleus(132,50);
        Nucleus projectile_h = get_default_nucleus(1,1);

        GlauberModelOLA_FM gm(projectile, target, range_t::FiniteRange);
        GlauberModelOLA_FM gmd(projectile_h, target, range_t::FiniteRange);

        range_integrated_type range_integrated;
        for(double b:{0.01,0.1, 0.2, 1.0})
            for(double e:{30.0,70.0,100.0,200.0})
            for(nucleon_t nn:{nucleon_t::pp, nucleon_t::pn, nucleon_t::np, nucleon_t::nn})
                {
                range_integrated.calculate(gm.projectile, gm.target, gm.beta());
                auto projectile_z = get_projectile(gm.z_integrated, nn);
                auto target_z = get_target(gm.z_integrated, nn);
                auto projectile_zr = get_projectile(range_integrated, nn);
                auto target_zr = get_target(range_integrated, nn);

                double p1 = ola_xyintegral_NN(projectile_z, target_z,b, gm.beta());
                double p2 = ola_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr, b);
                CHECK(p1 == approx(p2).R(rel_prec));
                // symmetry test
                double p3 = ola_xyintegral_NN(target_z, projectile_z, b, gm.beta());
                double p4 = ola_xyintegral_constrange_NN(target_z, projectile_z, target_zr, projectile_zr, b);
                CHECK(p3 == approx(p4).R(rel_prec));
                CHECK(p3 == approx(p1).R(rel_prec));
                CHECK(p4 == approx(p2).R(rel_prec));
            }
    }

TEST_CASE("mol4cfr_mol4cfr_constrange"){
        Nucleus projectile = get_default_nucleus(12,6);
        Nucleus target = get_default_nucleus(27,13);
        GlauberModel<MOL4C> gm(projectile, target, range_t::FiniteRange);
        Nucleus projectile_h = get_default_nucleus(1,1);
        GlauberModel<MOL4C> gmd(projectile_h, target, range_t::FiniteRange);
        range_integrated_type range_integrated;
        range_integrated.calculate(gm.projectile, gm.target, gm.beta());
        for(double b:{0.01,0.1, 0.2, 1.0})
            for(double e:{30.0,70.0,100.0,200.0}){                
                double sigmann = 0.5*gm.template sigmann<nucleon_t::pp>(e); // to be consistent profile function normalisation to 0.5
                auto projectile_z = get_projectile<nucleon_t::pp>(gm.z_integrated);
                auto target_z = get_target<nucleon_t::pp>(gm.z_integrated);
                auto projectile_zr = get_projectile<nucleon_t::pp>(range_integrated);
                auto target_zr = get_target<nucleon_t::pp>(range_integrated);

                double p1 = mol4c_xyintegral_NN<range_t::FiniteRange>(projectile_z, target_z, b, sigmann, gm.beta());
                double p2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z,
                                                            projectile_zr, target_zr,
                                                              b,sigmann);
               CHECK(p1 == approx(p2).epsilon(0.001));

                double p3 = mol4c_xyintegral_NN<range_t::FiniteRange>(projectile_z, target_z, b, sigmann, gm.beta());
                double p4 = mol4c_xyintegral_constrange_NN(projectile_z, target_z,
                                                            projectile_zr, target_zr,
                                                              b,sigmann);
               CHECK(p3 == approx(p4).epsilon(0.001));
                double p5 = mol4c_xyintegral_NN<range_t::FiniteRange>(projectile_z, target_z, b, sigmann, gm.beta());
                double p6 = mol4c_xyintegral_constrange_NN(projectile_z, target_z,
                                                            projectile_zr, target_zr,
                                                              b,sigmann);
               CHECK(p5 == approx(p6).epsilon(0.001));
            }
    }

TEST_CASE("ola_fmd"){
        const double rel_prec = 0.001;
        Nucleus projectile = get_default_nucleus(12,6);
        Nucleus target = get_default_nucleus(27,13);
        Nucleus projectile_h = get_default_nucleus(1,1);

        GlauberModel<OLA_FMD, FermiMotionD<NNCrossSectionFit>> gm1(projectile, target, range_t::FiniteRange);
        GlauberModel<OLA_FMD, FermiMotionD<NNCrossSectionFit>> gm2(target, projectile, range_t::FiniteRange);
        GlauberModel<OLA_FMD, FermiMotionD<NNCrossSectionFit>> gmd(projectile_h, target, range_t::FiniteRange);
        GlauberModel<OLA_FMD, FermiMotionD<NNCrossSectionFit>> gmd2(target, projectile_h, range_t::FiniteRange);        
        for(double b:{0.01,0.1, 0.2, 1.0})
            for(double e:{30.0,70.0,100.0,200.0}){
                double p1 = ola_fm_xyintegral<nucleon_t::pp>(gm1,b, e);
                double p2 = ola_fm_xyintegral_constrange<nucleon_t::pp>(gm1,b, e);
                CHECK(p1 == approx(p2).R(rel_prec));
                CHECK(p1 == approx(ola_fm_xyintegral<nucleon_t::pp>(gm2,b, e)).R(0.01));
                CHECK(p2 == approx(ola_fm_xyintegral_constrange<nucleon_t::pp>(gm2,b, e)).R(0.01));
            }
        for(double b:{0.01,0.1, 0.2, 1.0})
            for(double e:{30.0,70.0,100.0,200.0}){
                double p1 = ola_fm_xyintegral_dirac<nucleon_t::pp>(gmd,b, e);
                double p2 = ola_fm_xyintegral_dirac<nucleon_t::pp>(gmd2,b, e);
                CHECK(p1 == approx(p2).R(rel_prec));

            }
    }

TEST_CASE("mol4cfmd"){
        Nucleus projectile = get_default_nucleus(12,6);
        Nucleus target = get_default_nucleus(27,13);
        Nucleus projectile_h = get_default_nucleus(1,1);
        GlauberModel<MOL4C_FMD, FermiMotionD<NNCrossSectionFit>> gm1(projectile, target, range_t::FiniteRange);
        GlauberModel<MOL4C_FMD, FermiMotionD<NNCrossSectionFit>> gmd(projectile_h, target, range_t::FiniteRange);
        GlauberModel<MOL4C_FMD, FermiMotionD<NNCrossSectionFit>> gmd2(target, projectile_h, range_t::FiniteRange);

        for(double b:{0.01,0.1, 0.2, 1.0})
            for(double e:{30.0,70.0,100.0,200.0}){
                double p1 = mol4_fm_xyintegral<nucleon_t::pp,range_t::FiniteRange>(gm1,b, e);
                double p2 = mol4_fm_xyintegral_constrange<nucleon_t::pp>(gm1,b, e);
                CHECK(p1 == approx(p2).epsilon(0.001));
            }
        for(double b:{0.01,0.1, 0.2, 1.0})
            for(double e:{30.0,70.0,100.0,200.0}){
                double p1 = mol4_fm_xyintegral<nucleon_t::pp,range_t::FiniteRange>(gmd,b, e);
                double p2 = mol4_fm_xyintegral<nucleon_t::pp,range_t::FiniteRange>(gmd2,b, e);
                CHECK(p1 == approx(p2).epsilon(0.001));
            }
    }

TEST_CASE ("mol4c_vs_mol4cfmd"){
        Nucleus projectile = get_default_nucleus(12,6);
        Nucleus target = get_default_nucleus(27,13);
        GlauberModel<MOL4C_FMD, NNCrossSection_FermiMotionD> gm3(projectile, target,  range_t::FiniteRange);
        gm3.SetFECoefficient(0.01);
        GlauberModel<MOL4C_FMD> gm4(projectile, target,  range_t::FiniteRange);

        CHECK(gm3.beta() == gm4.beta());
        CHECK(gm3.beta() == gm4.beta());
        CHECK(gm3.fe_coefficient == approx(0.01,1e-4));

        double r1, r2;
        for(double b:{0.01,0.1, 0.2, 0.5, 1.0, 1.5})
            for(double e:{30.0,70.0,100.0,200.0,300.0}){

                auto projectile_z = get_projectile<nucleon_t::np>(gm4.z_integrated);
                auto target_z = get_target<nucleon_t::np>(gm4.z_integrated);
                auto projectile_zr = get_projectile<nucleon_t::np>(gm4.range_integrated);
                auto target_zr = get_target<nucleon_t::np>(gm4.range_integrated);
                double sigmann = 0.5*gm4.template sigmann<nucleon_t::np>(e);

                r1 = mol4_fm_xyintegral_constrange<nucleon_t::np>(gm3,b, e);
                r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
                CHECK(r1 == approx(r2).R(0.01));

                projectile_z = get_projectile<nucleon_t::nn>(gm4.z_integrated);
                target_z = get_target<nucleon_t::nn>(gm4.z_integrated);
                projectile_zr = get_projectile<nucleon_t::nn>(gm4.range_integrated);
                target_zr = get_target<nucleon_t::nn>(gm4.range_integrated);
                sigmann = 0.5*gm4.template sigmann<nucleon_t::nn>(e);

                r1 = mol4_fm_xyintegral_constrange<nucleon_t::nn>(gm3,b, e);
                r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
                CHECK(r1 == approx(r2).R(0.01));

                projectile_z = get_projectile<nucleon_t::pp>(gm4.z_integrated);
                target_z = get_target<nucleon_t::pp>(gm4.z_integrated);
                projectile_zr = get_projectile<nucleon_t::pp>(gm4.range_integrated);
                target_zr = get_target<nucleon_t::pp>(gm4.range_integrated);
                sigmann = 0.5*gm4.template sigmann<nucleon_t::pp>(e);
                r1 = mol4_fm_xyintegral_constrange<nucleon_t::pp>(gm3,b, e);
                r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
                CHECK(r1 == approx(r2).R(0.01));

                projectile_z = get_projectile<nucleon_t::pn>(gm4.z_integrated);
                target_z = get_target<nucleon_t::pn>(gm4.z_integrated);
                projectile_zr = get_projectile<nucleon_t::pn>(gm4.range_integrated);
                target_zr = get_target<nucleon_t::pn>(gm4.range_integrated);
                sigmann = 0.5*gm4.template sigmann<nucleon_t::pn>(e);
                r1 = mol4_fm_xyintegral_constrange<nucleon_t::pn>(gm3,b, e);
                r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
                CHECK(r1 == approx(r2).R(0.01));
            }

       gm3.SetFECoefficient(0.6);
       CHECK(gm3.fe_coefficient == approx(0.6,1e-4));

       // this is check at high energy, where FM should be negligible
       for(double b:{0.01,0.1, 0.2, 0.5, 1.0, 1.5,2.0,2.5})
           for(double e:{1500.,2000.}){
               auto projectile_z = get_projectile<nucleon_t::np>(gm4.z_integrated);
               auto target_z = get_target<nucleon_t::np>(gm4.z_integrated);
               auto projectile_zr = get_projectile<nucleon_t::np>(gm4.range_integrated);
               auto target_zr = get_target<nucleon_t::np>(gm4.range_integrated);
               double sigmann = 0.5*gm4.template sigmann<nucleon_t::np>(e);

               r1 = mol4_fm_xyintegral_constrange<nucleon_t::np>(gm3,b, e);
               r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
               CHECK(r1 == approx(r2).R(0.01));

               projectile_z = get_projectile<nucleon_t::nn>(gm4.z_integrated);
               target_z = get_target<nucleon_t::nn>(gm4.z_integrated);
               projectile_zr = get_projectile<nucleon_t::nn>(gm4.range_integrated);
               target_zr = get_target<nucleon_t::nn>(gm4.range_integrated);
               sigmann = 0.5*gm4.template sigmann<nucleon_t::nn>(e);

               r1 = mol4_fm_xyintegral_constrange<nucleon_t::nn>(gm3,b, e);
               r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
               CHECK(r1 == approx(r2).R(0.01));

               projectile_z = get_projectile<nucleon_t::pp>(gm4.z_integrated);
               target_z = get_target<nucleon_t::pp>(gm4.z_integrated);
               projectile_zr = get_projectile<nucleon_t::pp>(gm4.range_integrated);
               target_zr = get_target<nucleon_t::pp>(gm4.range_integrated);
               sigmann = 0.5*gm4.template sigmann<nucleon_t::pp>(e);
               r1 = mol4_fm_xyintegral_constrange<nucleon_t::pp>(gm3,b, e);
               r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
               CHECK(r1 == approx(r2).R(0.01));

               projectile_z = get_projectile<nucleon_t::pn>(gm4.z_integrated);
               target_z = get_target<nucleon_t::pn>(gm4.z_integrated);
               projectile_zr = get_projectile<nucleon_t::pn>(gm4.range_integrated);
               target_zr = get_target<nucleon_t::pn>(gm4.range_integrated);
               sigmann = 0.5*gm4.template sigmann<nucleon_t::pn>(e);
               r1 = mol4_fm_xyintegral_constrange<nucleon_t::pn>(gm3,b, e);
               r2 = mol4c_xyintegral_constrange_NN(projectile_z, target_z, projectile_zr, target_zr,b,sigmann);
               CHECK(r1 == approx(r2).R(0.01));
           }
    }

TEST_CASE("mol_constrange"){
    Nucleus projectile = get_default_nucleus(12,6);
    Nucleus target = get_default_nucleus(27,13);
    Nucleus projectile_h = get_default_nucleus(1,1);

    GlauberModel<MOL> gm1(projectile, target, range_t::FiniteRange);
    GlauberModel<MOL> gmd(projectile_h, target, range_t::FiniteRange);
    for(double b:{0.01,0.1, 0.2, 1.0})
        for(double E:{30.0,70.0,100.0,200.0}){

            double p1 = mol_xyintegral<nucleonmol_t::projectile_p,range_t::FiniteRange>(gm1,b,E);
            double p2 = mol_xyintegral_constrange<nucleonmol_t::projectile_p,range_t::FiniteRange>(gm1,b,E);
            CHECK(p1 == approx(p2).epsilon(0.001));

            double p3 = mol_xyintegral<nucleonmol_t::projectile_n,range_t::FiniteRange>(gm1,b,E);
            double p4 = mol_xyintegral_constrange<nucleonmol_t::projectile_n,range_t::FiniteRange>(gm1,b,E);
            CHECK(p3 == approx(p4).epsilon(0.001));

            double p5 = mol_xyintegral<nucleonmol_t::target_p,range_t::FiniteRange>(gm1,b,E);
            double p6 = mol_xyintegral_constrange<nucleonmol_t::target_p,range_t::FiniteRange>(gm1,b,E);
            CHECK(p5 == approx(p6).epsilon(0.001));

            double p7 = mol_xyintegral<nucleonmol_t::target_n,range_t::FiniteRange>(gm1,b,E);
            double p8 = mol_xyintegral_constrange<nucleonmol_t::target_n,range_t::FiniteRange>(gm1,b,E);
            CHECK(p7 == approx(p8).epsilon(0.001));
        }

    /*
    for(double b:{0.01,0.1, 0.2, 1.0})
        for(double E:{30.0,70.0,100.0,200.0}){
            double p1 = mol_xyintegral<nucleonmol_t::projectile_p,range_t::FiniteRange>(gmd,b,E);
            double p2 = mol_xyintegral_constrange<nucleonmol_t::projectile_p,range_t::FiniteRange>(gmd,b,E);
            TEST_CHECK(p1 == approx(p2).epsilon(0.001));*/
            /*
            double p3 = mol_xyintegral<nucleon_t::target_n,range_t::FiniteRange>(gmd,b,E);
            double p4 = mol_xyintegral_constrange<nucleon_t::target_n,range_t::FiniteRange>(gmd,b,E);
            TEST_CHECK(p3 == approx(p4).epsilon(0.001));
            double p5 = mol_xyintegral<nucleon_t::target_n,range_t::FiniteRange>(gmd,b,E);
            double p6 = mol_xyintegral_constrange<nucleon_t::target_n,range_t::FiniteRange>(gmd,b,E);
            TEST_CHECK(p5 == approx(p6).epsilon(0.001));*/
        //}
}

TEST_CASE("mol_vs_molfmd"){
        Nucleus projectile = get_default_nucleus(12,6);
        Nucleus target = get_default_nucleus(27,13);
        GlauberModel<MOL_FMD, NNCrossSection_FermiMotionD> gm3(projectile, target,  range_t::FiniteRange);
        gm3.SetFECoefficient(0.01);
        GlauberModel<MOL> gm4(projectile, target,  range_t::FiniteRange);
        for(double b:{0.01,0.1, 0.2, 0.5, 1.0, 1.5})
            for(double e:{30.0,70.0,100.0,200.0,300.0}){
                double r1,r2;
                r1 = mol_xyintegral_constrange_fm<nucleonmol_t::projectile_p, range_t::FiniteRange>(gm3,b, e);
                r2 = mol_xyintegral_constrange<nucleonmol_t::projectile_p, range_t::FiniteRange>(gm4,b,e);
                CHECK(r1 == approx(r2).R(0.01));
                r1 = mol_xyintegral_constrange_fm<nucleonmol_t::projectile_n, range_t::FiniteRange>(gm3,b, e);
                r2 = mol_xyintegral_constrange<nucleonmol_t::projectile_n, range_t::FiniteRange>(gm4,b,e);
                CHECK(r1 == approx(r2).R(0.01));
                r1 = mol_xyintegral_constrange_fm<nucleonmol_t::target_n, range_t::FiniteRange>(gm3,b, e);
                r2 = mol_xyintegral_constrange<nucleonmol_t::target_n, range_t::FiniteRange>(gm4,b,e);
                CHECK(r1 == approx(r2).R(0.01));
                r1 = mol_xyintegral_constrange_fm<nucleonmol_t::target_p, range_t::FiniteRange>(gm3,b, e);
                r2 = mol_xyintegral_constrange<nucleonmol_t::target_p, range_t::FiniteRange>(gm4,b,e);
                CHECK(r1 == approx(r2).R(0.01));
            }

       gm3.SetFECoefficient(0.6);
       for(double b:{0.01,0.1, 0.2, 0.5, 1.0, 1.5,2.0,2.5})
           for(double e:{1500.,2000.}){
               double r1, r2;
               r1 = mol_xyintegral_constrange_fm<nucleonmol_t::projectile_p, range_t::FiniteRange>(gm3,b, e);
               r2 = mol_xyintegral_constrange<nucleonmol_t::projectile_p, range_t::FiniteRange>(gm4,b,e);
               CHECK(r1 == approx(r2).R(0.01));
               r1 = mol_xyintegral_constrange_fm<nucleonmol_t::projectile_n, range_t::FiniteRange>(gm3,b, e);
               r2 = mol_xyintegral_constrange<nucleonmol_t::projectile_n, range_t::FiniteRange>(gm4,b,e);
               CHECK(r1 == approx(r2).R(0.01));
               r1 = mol_xyintegral_constrange_fm<nucleonmol_t::target_n, range_t::FiniteRange>(gm3,b, e);
               r2 = mol_xyintegral_constrange<nucleonmol_t::target_n, range_t::FiniteRange>(gm4,b,e);
               CHECK(r1 == approx(r2).R(0.01));
               r1 = mol_xyintegral_constrange_fm<nucleonmol_t::target_p, range_t::FiniteRange>(gm3,b, e);
               r2 = mol_xyintegral_constrange<nucleonmol_t::target_p, range_t::FiniteRange>(gm4,b,e);
               CHECK(r1 == approx(r2).R(0.01));
           }
    }


