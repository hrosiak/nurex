#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include "nurex/NNCrossSection.h"
#include "nurex/nurex.h"

using namespace nurex;

    TEST_CASE("Fermi Motion Correction"){
        FermiMotion<NNCrossSectionFit> sigmaNN(90,90);
        
        CHECK(sigmaNN.np(0.1)==approx(3895,5));
        CHECK(sigmaNN.pp(0.1)==approx(377,5));
        
        CHECK(sigmaNN.np(1)==approx(3735,5));
        CHECK(sigmaNN.pp(1)==approx(369,5));
        
        CHECK(sigmaNN.np(10)==approx(2472,5));
        CHECK(sigmaNN.pp(10)==approx(296,5));
        
        CHECK(sigmaNN.np(30)==approx(1039,3));
        CHECK(sigmaNN.pp(30)==approx(177,3));
        
        CHECK(sigmaNN.np(50)==approx(482,2));
        CHECK(sigmaNN.pp(50)==approx(108.3,2));
        
        CHECK(sigmaNN.np(100)==approx(125.8,1));
        CHECK(sigmaNN.pp(100)==approx(42.6,1));
        
        CHECK(sigmaNN.np(300)==approx(37.5,1));
        CHECK(sigmaNN.pp(300)==approx(24.5,1));
        
        CHECK(sigmaNN.np(650)==approx(36.5,1));
        CHECK(sigmaNN.pp(650)==approx(40.0,1));
        
        NNCrossSectionFit sigma_nn;
        CHECK(sigmaNN.np(2000)==approx(sigma_nn.np(2000),1));
        CHECK(sigmaNN.pp(2000)==approx(sigma_nn.pp(2000),1));  
    }

    TEST_CASE("zero correction"){
        FermiMotion<NNCrossSectionFit> sigmaNN(0,0);
        NNCrossSectionFit sigma_nn;
        for(double e:{1,10,50,100,500,600}){
            CHECK(sigmaNN.np(e)==approx(sigma_nn.np(e),1));
            CHECK(sigmaNN.pp(e)==approx(sigma_nn.pp(e),1));
        }
    }

    TEST_CASE("negative momentum"){
        NNCrossSectionFit sigma_nn;

        for(double p=0.01; p<p_from_T(2000);p*=1.2){
            CHECK_MESSAGE(sigma_nn.np(T_from_p(p))>0,"E = "<<T_from_p(p));
            CHECK(sigma_nn.np(T_from_p(p)) == sigma_nn.np(T_from_p(-p)));
            CHECK(sigma_nn.pp(T_from_p(p)) == sigma_nn.pp(T_from_p(-p)));
        }
    }

    TEST_CASE("copy test"){
        FermiMotion<NNCrossSectionFit> sigma_nn(90,18);
        FermiMotion<NNCrossSectionFit> s2 = sigma_nn;
        
        CHECK(sigma_nn.np(30)==approx(s2.np(30),0.001));
        CHECK(sigma_nn.np(100)==approx(s2.np(100),0.001));
        CHECK(sigma_nn.pp(30)==approx(s2.pp(30),0.001));
        CHECK(sigma_nn.pp(100)==approx(s2.pp(100),0.001));

        sigma_nn.SetMomentum(90,90);
        CHECK(sigma_nn.np(30)!=approx(s2.np(30),0.001));
        CHECK(sigma_nn.np(100)!=approx(s2.np(100),0.001));
        CHECK(sigma_nn.pp(30)!=approx(s2.pp(30),0.001));
        CHECK(sigma_nn.pp(100)!=approx(s2.pp(100),0.001));
    }
    
     TEST_CASE("non constant FermiMotion"){
         FermiMotion<NNCrossSectionFit> csnn(90,90);
         FermiMotionD<NNCrossSectionFit> snn;
         NNCrossSectionFit s0;
         Nucleus p= get_default_nucleus(12,6);
        
         for(double e : {10,20,30,50,100,200,600}){
             CHECK(csnn.pp(e)==approx(snn.pp(e,90,90),1));
             CHECK(csnn.np(e)==approx(snn.np(e,90,90),1));
             CHECK(s0.pp(e)==approx(snn.pp(e,0,0),0.1));
             CHECK(s0.np(e)==approx(snn.np(e,0,0),0.1));
         }
    }