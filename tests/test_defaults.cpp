#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include <vector>
#include "nurex/Nucleus.h"
#include "nurex/DefaultDensities.h"
using namespace nurex;

#define TEST_CHECK CHECK

TEST_CASE("default nucleus"){
        Nucleus p = get_default_nucleus(12,6);
        TEST_CHECK(p.GetDensityProton().Norm() == 6);
        TEST_CHECK(p.A() == 12);
        TEST_CHECK(p.RrmsProton()==approx(2.322).epsilon(0.05));
        TEST_CHECK(p.RrmsNeutron()==approx(2.322).epsilon(0.1));;
        
        p = get_default_nucleus(60,28);
        TEST_CHECK(!(!p));
	    TEST_CHECK(p.A() == 60);
        TEST_CHECK(p.GetDensityProton().Norm() == 28);

        p = get_default_nucleus(6,1);
        TEST_CHECK(!p);

        Nucleus p2 = get_default_nucleus("12c");
        TEST_CHECK(p2.GetDensityProton().Norm() == 6);
        TEST_CHECK(p2.A() == 12);
        TEST_CHECK(p2.RrmsProton()==approx(2.322).epsilon(0.05));
        TEST_CHECK(p2.RrmsNeutron()==approx(2.322).epsilon(0.1));
    }
    
TEST_CASE("special"){
        Nucleus p = get_default_nucleus("1H");
        Nucleus p2 = get_default_nucleus(1,1);
        TEST_CHECK(p.A() == 1);
        TEST_CHECK(p.Z() == 1);
        TEST_CHECK(p.N() == 0);
        TEST_CHECK(p.GetDensityProton().Norm() == 1);
        TEST_CHECK(p.GetDensityNeutron().Norm() == 0);
        TEST_CHECK(p.RrmsNeutron()==0.0);
        TEST_CHECK(p.RrmsProton()<0.001);
        TEST_CHECK(p2.A() == 1);
        TEST_CHECK(p2.Z() == 1);
        TEST_CHECK(p2.N() == 0);
        TEST_CHECK(p2.GetDensityProton().Norm() == 1);
        TEST_CHECK(p2.GetDensityNeutron().Norm() == 0);
    }    

