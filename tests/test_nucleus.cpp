#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include <vector>
#include "nurex/Nucleus.h"
#include "nurex/DefaultDensities.h"

using namespace nurex;

    TEST_CASE("Nucleus Check"){
    
        auto df0 = make_density<DensityFermi>(4.2413494,0.50732378);
        auto df1 = make_density<DensityFermi>(4.0,0.5);

        Nucleus projectile(62,28);
        CHECK(!projectile);
        CHECK(!projectile.GetDensityNeutron());
        CHECK(!projectile.GetDensityProton());
        
        Nucleus projectile2(62,28, df0, df1);
        CHECK_FALSE(!projectile2);
        CHECK(projectile2.GetDensityProton().Norm()==approx(28).epsilon(1e-3));
        CHECK(projectile2.GetDensityNeutron().Norm()==approx(34).epsilon(1e-3));
        CHECK(projectile2.DensityProton(0)==approx(28*0.00274114).epsilon(1e-6));
        
        Nucleus proj2(62,28,DensityFermi(4.2413494,0.50732378), DensityFermi(4.2413494,0.50732378));
        CHECK(proj2.GetDensityProton().Norm() == 28);
        CHECK(proj2.GetDensityNeutron().Norm() == 34);
        
    }
    
    TEST_CASE("Nucleus(references)"){
        auto df1 = make_density<DensityFermi>(4.24,0.507);
        auto df2 = make_density<DensityFermi>(4.24,0.507);
        Nucleus proj(62,28,df1,df2);
        CHECK(proj.GetDensityProton().Norm() == 28);
        CHECK(proj.GetDensityNeutron().Norm() == 34);
    }
    
    TEST_CASE("Nucleus move constructor"){
        Nucleus proj(get_default_nucleus(12,6));
        CHECK(proj.GetDensityProton().Norm() == 6);
        CHECK(proj.GetDensityNeutron().Norm() == 6);
    }
    
    TEST_CASE("density copy when initialized"){
        auto df1 = make_density<DensityFermi>(4.24,0.507);
        auto df2 = make_density<DensityFermi>(4.24,0.507);
        Nucleus proj(62,28,df1,df2);
        double r1 = proj.RrmsProton();
        df1.SetParameters(6,0.7);
        double r2 = proj.RrmsProton();
        CHECK(r2 == r1);
        CHECK(Rrms(df1)>r2);
        CHECK(df1 != proj.GetDensityProton());
    }
    
    TEST_CASE("not properly initialized"){
        Nucleus p1(62,28);
        Nucleus p2 = get_default_nucleus(6,1);
        Nucleus p3;
        
        CHECK_FALSE(p1);
        CHECK(!p2);
        CHECK(!p3);
    }
    
    TEST_CASE("Nucleus without neutron"){
        Nucleus proj(1,1,make_density<DensityFermi>(4.24,0.507),make_density<DensityZero>());
        CHECK(proj.N()==0);
        CHECK(proj.RrmsNeutron()==0.0);
    }
    
    TEST_CASE("Nucleus string init"){
        Nucleus p("12c");
        CHECK(p.GetDensityProton().Norm() == 6);
        CHECK(p.GetDensityNeutron().Norm() == 6);
        CHECK(p.RrmsProton()==approx(2.322).epsilon(0.001));
        CHECK(p.RrmsNeutron()==approx(2.322).epsilon(0.1));
        CHECK_THROWS(Nucleus p("7c"));
    }
    
   TEST_CASE("Nucleus copy"){
        Nucleus p1(62,28,make_density<DensityFermi>(4.24,0.507),make_density<DensityFermi>(4.24,0.507));
        Nucleus p2 = p1;
        CHECK(p2.GetDensityProton().Norm() == 28);
        CHECK(p2.GetDensityNeutron().Norm() == 34);
        CHECK(p2.A() == p1.A());
        CHECK(p2.Z() == p1.Z());
        
        CHECK( p2.GetDensityProton().Rrms() == approx(p1.GetDensityProton().Rrms()).epsilon(1e-6));
        CHECK( p2.GetDensityNeutron().Rrms() == approx(p1.GetDensityNeutron().Rrms()).epsilon(1e6));
    }
    
    TEST_CASE("Nucleus move"){
        std::vector<Nucleus> v{Nucleus(62,28,make_density<DensityFermi>(4.24,0.507),make_density<DensityFermi>(4.24,0.507)),Nucleus(62,28,make_density<DensityFermi>(4.24,0.507),make_density<DensityFermi>(4.24,0.507))};
        CHECK(v.size()==2);
        CHECK(v[0].GetDensityProton().Norm() == 28);
        CHECK(v[1].GetDensityProton().Norm() == 28);
    }
    
    TEST_CASE("equality"){
        Nucleus p1(62,28);
        Nucleus p2 = get_default_nucleus(12,6);
        Nucleus p3;
        CHECK(p1!=p2);
        CHECK(p1!=p3);
        CHECK(p2!=p3);

        auto def1 = get_default_nucleus(12,6);
        CHECK(def1==p2);

        auto df1 = make_density<DensityFermi>(4.24,0.507);
        auto df2 = make_density<DensityFermi>(4.24,0.507);
        auto df3 = make_density<DensityFermi>(4.24,0.57);
        CHECK(df1==df2);
        CHECK(df1!=df3);
        Nucleus p4(62,28,df1,df2);
        Nucleus p5(62,28,df1,df2);
        Nucleus p6(62,29,df1,df2);
        Nucleus p7(62,28,df1,df3);
        CHECK(p4==p5);
        CHECK(p4!=p6);
        CHECK(p4!=p7);
    }
    
    TEST_CASE("compounds"){
        Nucleus c = get_default_nucleus(12,6);
        Nucleus h = get_default_nucleus(1,1);       
        Nucleus o = get_default_nucleus(16,8);       
        Compound ch2 = Compound({{c,1},{h,2}});
        
        CHECK(ch2.ncomponents()==2);
        CHECK(ch2.get_nucleus(0)==c);
        CHECK(ch2.get_nucleus(1)==h);
        CHECK(ch2.molar_fraction(0)==1);
        CHECK(ch2.molar_fraction(1)==2);
        
        Compound water = Compound({
                {{"1h"},2},
                {{"16o"},1}
        });
        CHECK(water.ncomponents()==2);
        CHECK(water.get_nucleus(0)==h);
        CHECK(water.get_nucleus(1)==o);
        
        auto my_o = Nucleus(16,8,make_density<DensityFermi>(1.5,0.7),make_density<DensityFermi>(1.5,0.7));
        Compound water2= Compound({
                {h,2},
                {Nucleus(16,8,make_density<DensityFermi>(1.5,0.7),make_density<DensityFermi>(1.5,0.7)),1}
        });
        CHECK(water2.ncomponents()==2);
        CHECK(water2.get_nucleus(0)==h);
        CHECK(water2.get_nucleus(1)==my_o);

        }
        

