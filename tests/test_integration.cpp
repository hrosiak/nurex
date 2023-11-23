#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include "nurex/nurex.h"

using namespace nurex;

    TEST_CASE("Z Integration test"){
    DensityType df0 = DensityFermi(4.2413494,0.50732378);
    df0.Normalize(62);
    //Interpolator *zint = dynamic_cast<Interpolator*>(ZIntegrate(&df0));
    Functional zint = ZIntegrate(df0);
    CHECK(zint.eval(0)==approx(1.441733109).epsilon(1e-3));
    CHECK(zint.eval(1.0023322)==approx(1.3981238477).epsilon(1e-3));
    CHECK(zint.eval(3.66)==approx(0.68365512142).epsilon(1e-3));
    CHECK(zint.eval(5.49)==approx(0.059374015325).epsilon(1e-3));
    CHECK(zint.min() == approx(0.0).epsilon(1e-6));
    CHECK(zint.max() == approx(max_rrms_factor*df0.Rrms()).epsilon(1e-6));
    CHECK_FALSE(is_ftype<ConstantFunction>(zint));
    CHECK_FALSE(zint.is_type<ConstantFunction>());

    DensityType dz = DensityZero();
    Functional zint2 = ZIntegrate(dz);

    CHECK(is_ftype<ConstantFunction>(zint2));
    CHECK(zint2.is_type<ConstantFunction>());

    CHECK(zint2.eval(0.0)==0.0);
    CHECK(zint2.eval(0.0001)==0.0);

    DensityType dd = DensityDirac();
    Functional zint3 = ZIntegrate(dd);

    CHECK(is_ftype<DiracFunction>(zint3));
    CHECK_FALSE(is_ftype<ConstantFunction>(zint3));
    CHECK(zint3.is_type<DiracFunction>());
    CHECK_FALSE(zint3.is_type<ConstantFunction>());
    }
    
    TEST_CASE("2d norm"){
        double i;
        DensityType df = make_density<DensityFermi>(1,0.5);
        auto zf1 = ZIntegrate(df);
        auto f1 = [&](double x, double y){
            return zf1.eval(sqrt(x*x+y*y));
        };
        double rmax = Rrms(df)*4;
        GaussLegendreIntegration2D<32> int2d;
        i = int2d(f1,-rmax,rmax,-rmax,rmax);
        CHECK(i==approx(1.0).epsilon(0.0005));
    }
    
    TEST_CASE("2d norm gaussians"){
        double i;
        DensityType df = make_density<DensityGaussian>(1);
        auto zf1 = ZIntegrate(df);
        auto f1 = [&](double x, double y){
            return zf1.eval(sqrt(x*x+y*y));
        };
        double rmax = Rrms(df)*4;
        GaussLegendreIntegration2D<32> int2d;
        i = int2d(f1,-rmax,rmax,-rmax,rmax);
        CHECK(i==approx(1.0).epsilon(0.0001));
    }
    
    TEST_CASE("Finite Range Normalization check"){
        auto f = [](double x, double y){
            double r = sqrt(x*x + y*y);
            return finite_range(r,0.39);
        };
        
        auto f2 = [](double x, double y){
            double r = sqrt(x*x + y*y);
            return finite_range(r,0.039);
        };
        
        auto f3 = [](double x, double y){
            double r = sqrt(x*x + y*y);
            return finite_range(r,0.99);
        };
        GaussLegendreIntegration2D<32> int2d;
        double i = int2d(f,-3,3,-5,5);
        CHECK(i == approx(1.0).epsilon(1e-2));
        i = int2d(f2,-0.2,0.21,-0.21,0.21);
        CHECK(i == approx(1.0).epsilon(1e-2));
        i = int2d(f3,-5,5,-5,5);
        CHECK(i == approx(1.0).epsilon(1e-2));
    }