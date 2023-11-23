#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include "nurex/Density.h"
#include "nurex/Nucleus.h"
using namespace nurex;
#define TEST_CHECK CHECK

TEST_CASE("normalisation"){
        DensityFermi df1(4.2413494,0.50732378);
        DensityFermi df2(4.2413494,0.50732378);
        TEST_CHECK(df1.Norm()==1);
        df1.Normalize(26);
        normalize(df2,26);
        TEST_CHECK(df1.Norm()==26);
        TEST_CHECK(df2.Norm()==26);
        df1.Normalize();
        normalize(df2);
        TEST_CHECK(df1.Norm()==26);
        TEST_CHECK(df2.Norm()==26);

        normalize(df1,34);
        TEST_CHECK(df1.Norm()==34);


        // constructor with normalization
        DensityFermi df4(4.2413494,0.50732378,0, 26);
        TEST_CHECK(df4.Norm()==26);
        df4.Normalize();
        TEST_CHECK(df4.Norm()==26);
        df4.Normalize(10);
        TEST_CHECK(df4.Norm()==10);

        DensityHO df5(1.6,1.6,26);
        TEST_CHECK(df5.Norm()==26);

    }
    
    
TEST_CASE("fermi"){
        auto df1m = make_density<DensityFermi>(4.2413494,0.50732378);
        auto df2m = make_density<DensityFermi>(4.2413494,0.50732378,0.1);
        auto df3 = make_density<DensityFermi>(4.2413494,0.50732378,0.1, 6);
        DensityFermi df1(4.2413494,0.50732378);
        DensityFermi df2(4.2413494,0.50732378,0.1);
        DensityFermi df4(4.2413494,0.50732378,0.001);

        TEST_CHECK(df1 == df1m);
        TEST_CHECK(df2 == df2m);
        TEST_CHECK(df3 != df2m);
        df2m.Normalize(6);
        TEST_CHECK(df3 == df2m);

        TEST_CHECK(df1.GetParameter(0)==4.2413494);
        TEST_CHECK(df1.GetParameter(1)==0.50732378);
        TEST_CHECK(df1.GetParameter(2)==0.0);
        TEST_CHECK(df1.nparameters==3);

        TEST_CHECK(df2.GetParameter(0)==4.2413494);
        TEST_CHECK(df2.GetParameter(1)==0.50732378);
        TEST_CHECK(df2.GetParameter(2)==0.1);
        TEST_CHECK(df2.nparameters==3);

        TEST_CHECK(df1(0)==approx(0.00274114,1e-6));
        TEST_CHECK(df1(0)!=approx(df2(0),1e-8));
        TEST_CHECK(df1(0.1)!=approx(df2(0.1),1e-8));
        TEST_CHECK(df1(0)==approx(df4(0),1e-4));
        TEST_CHECK(df1(0.1)==approx(df4(0.1),1e-4));

        TEST_CHECK(df1(2)==approx(0.00270912,1e-6));
        TEST_CHECK(df1(4)==approx(0.00169097,1e-6));
        TEST_CHECK(df1(5)==approx(0.00050206,1e-6));
        TEST_CHECK(df1(5)==approx(0.00050206,1e-6));
        TEST_CHECK(Rrms(df1)==approx(3.78810001,1e-6));

        df1.Normalize(100);
        TEST_CHECK(df1.GetParameter(0)==4.2413494);
        TEST_CHECK(df1.GetParameter(1)==0.50732378);
        TEST_CHECK(df1.GetParameter(2)==0.0);
        TEST_CHECK(df1(0)==approx(100*0.00274114,1e-6));
        TEST_CHECK(df1(2)==approx(100*0.00270912,1e-6));
        TEST_CHECK(df1(4)==approx(100*0.00169097,1e-6));
        TEST_CHECK(df1(5)==approx(100*0.00050206,1e-6));
        TEST_CHECK(Rrms(df1)==approx(3.78810001,1e-6));
    }

TEST_CASE("ho"){
        DensityHO df2(1.5484, 1.6038);

        TEST_CHECK(df2.GetParameter(0)==1.5484);
        TEST_CHECK(df2.GetParameter(1)==1.6038);
        TEST_CHECK(df2.GetParameter(2)==0.0);
        TEST_CHECK(df2.nparameters==2);

        TEST_CHECK(Rrms(df2)==approx(2.3,0.01));
    }

TEST_CASE("zero"){
        DensityZero dz;
        TEST_CHECK(dz.Density(0)==0.0);
        TEST_CHECK(Rrms(dz)==0.0);
        TEST_CHECK(dz.GetParameter(0)==0.0);
        TEST_CHECK(dz.GetParameter(1)==0.0);
        TEST_CHECK(dz.nparameters==0);

        auto dz2 = make_density<DensityZero>();
        TEST_CHECK(Rrms(dz2)==0.0);
        TEST_CHECK(dz2.Norm()==0.0);
    }
    
TEST_CASE("dirac"){
        DensityDirac dd;

        TEST_CHECK(dd.GetParameter(0)==0.0);
        TEST_CHECK(dd.GetParameter(1)==0.0);
        TEST_CHECK(dd.nparameters==0);

        TEST_CHECK(dd.Density(0.0)==0.0);
        TEST_CHECK(dd.Density(0.0001)==0.0);
        TEST_CHECK(dd.Norm()==1.0);
        dd.Normalize(2.0);
        TEST_CHECK(dd.Norm()==2.0);
        TEST_CHECK(dd.Density(0.0)==0.0);
    }
    
TEST_CASE("gauss"){
        DensityGaussian dg(1.0);

        TEST_CHECK(dg.GetParameter(0)==1.0);
        TEST_CHECK(dg.GetParameter(1)==0.0);
        TEST_CHECK(dg.nparameters==1);

        TEST_CHECK(dg.Norm()==approx(1.0).epsilon(1e-5));
        TEST_CHECK(dg.Density(0.0)==approx(0.0634936).epsilon(1e-5));
        TEST_CHECK(dg.GetRho0()==approx(0.0634936).epsilon(1e-5));
        dg.Normalize(2.0);
        TEST_CHECK(dg.Norm()==approx(2.0).epsilon(1e-5));
        TEST_CHECK(dg.Density(0.0)==approx(2*0.0634936).epsilon(1e-5));
        TEST_CHECK(dg.GetRho0()==approx(2*0.0634936).epsilon(1e-5));

    }
    
TEST_CASE("file"){
        #ifndef NO_FILESYSTEM
        bool r = false;
        try{
            density_from_file("nonexistingfile.txt");
            }
        catch(...){
            r = true;
            }
        TEST_CHECK(r);
        #endif
}
    
TEST_CASE("make"){
        auto df = make_density<DensityFermi>(4.2413494,0.50732378);
        auto df2 = make_density<DensityHO>(1.5484, 1.6038);
        TEST_CHECK(Rrms(df) == approx(3.78810001,1e-3));
        TEST_CHECK(df.type==density_type::fermi);
        TEST_CHECK(Rrms(df2) == approx(2.3,1e-3));
        TEST_CHECK(df2.type==density_type::ho);
    }
    
TEST_CASE("rrms"){
        auto df1 = make_density<DensityFermi>(4.09511,0.555);
        auto df2 = make_density<DensityFermi>(4.09713,0.555);

        double r1;
        double r2;
        r1 = Rrms(df1);
        df1.Normalize(26);
        r2 = Rrms(df1);
        TEST_CHECK(r1==approx(r2).epsilon(1e-5));

        df2.Normalize(26);
        r2 = Rrms(df2);
        TEST_CHECK(r1<r2);
        TEST_CHECK(df1.Density(0)>df2.Density(0));

        auto df3 = make_density<DensityFermi>(4.09713,0.555);
        //df3->Normalize();
        r1 = Rrms(df3);
        r2 = Rrms(df3);
        TEST_CHECK(r1==approx(r2).epsilon(1e-5));
        df3.SetParameters(df3.GetRadiusParameter()+0.1,df3.GetDiffusenessParameter());
        //df3->Normalize();
        double r3  = Rrms(df3);
        r2 = Rrms(df3);
        TEST_CHECK(r1<r3);
        TEST_CHECK(r1<r2);
        TEST_CHECK(r3==approx(r2).epsilon(1e-5));

    }

TEST_CASE("clone"){
        auto df1 = make_density<DensityFermi>(4.09511,0.555,64);
        auto df2 = df1;
        TEST_CHECK(Rrms(df1)==Rrms(df2));
        TEST_CHECK(df1.Norm()==df2.Norm());
        TEST_CHECK(df1.Density(0.0)==df2.Density(0.0));
        TEST_CHECK(df1.Density(0.5)==df2.Density(0.5));

        DensityType db1 = df1;
        DensityType db2 = db1;
        db1.Normalize(2.0);
        TEST_CHECK(0.5*db1.Density(0.0)==approx(db2.Density(0.0)).R(1e-3));
        TEST_CHECK(0.5*db1.Density(0.5)==approx(db2.Density(0.5)).R(1e-3));

      }
      
TEST_CASE("equality"){
        auto df1 = make_density<DensityFermi>(4.09511,0.555);
        auto df2 = df1;
        TEST_CHECK(df1== df2);
        auto df3 = make_density<DensityFermi>(4.09511,0.555);
        TEST_CHECK(df1 == df3);
        df3.Normalize(2.0);
        TEST_CHECK(df1 != df3);

        auto df4 = make_density<DensityFermi>(4.09511,0.555);
        auto df5 = make_density<DensityFermi>(4.09512,0.555);
        df4.Normalize(2.0);
        TEST_CHECK(df4 == df3);
        TEST_CHECK(df4 != df1);
        TEST_CHECK(df4 != df5);
        df5.Normalize(2.0);
        TEST_CHECK(df4 != df5);

        auto dh1 = make_density<DensityHO>(4.09511,0.555);
        auto dh2 = make_density<DensityHO>(4.09511,0.555);
        auto dh3 = make_density<DensityHO>(4.09511,0.556);
        TEST_CHECK(dh1 != df1);
        TEST_CHECK(dh1 != df3);
        TEST_CHECK(dh1 != df4);
        TEST_CHECK(dh1 == dh2);
        TEST_CHECK(dh1 != dh3);
    };
    
TEST_CASE("table"){
        auto df1 = make_density<DensityFermi>(4.09511,0.555,64);
        auto b = linspace_vector(0,4*Rrms(df1),100);

        std::vector<double> d;
        d.reserve(100);

        for(double x:b){
          d.push_back(df1.Density(x));
        }

        DensityTable dt(b,d);
        DensityTable dt2 = dt;
        TEST_CHECK(dt.Norm() == approx(1.0).R(1e-6));

        auto b1 = linspace_vector(0,2.0*Rrms(df1),10);
        auto b2 = linspace_vector(2.0,4.1*Rrms(df1),10);
        for(double x:b1){
          TEST_CHECK(dt.Density(x) == approx(df1(x)).epsilon(1e-7));
        }
        for(double x:b2){
          TEST_CHECK(dt.Density(x) == approx(df1(x)).epsilon(1e-7));
        }

        dt.Normalize(10.0);
        dt2.Normalize(2.0);
        for(double x:b1){
          TEST_CHECK(dt.Density(x) == approx(10.0*df1(x)).epsilon(1e-6));
          TEST_CHECK(dt2.Density(x) == approx(2.0*df1(x)).epsilon(1e-6));
        }

        Nucleus n1(12,6,dt,DensityZero());
        TEST_CHECK(!(!n1));

      }
      
TEST_CASE("string"){
        auto df1 = make_density<DensityFermi>(4.09511,0.555,64);
        TEST_CHECK(density_type_to_string(df1.type).compare("fermi")==0);
    }
    
TEST_CASE("base"){
    DensityFermi df1(4.2413494,0.50732378);
    DensityFermi df2(4.2413494,0.50732378);
    DensityType db1 = df1;
    DensityType db2(df1);

    TEST_CHECK(Rrms(db1)==approx(3.78810001,1e-6));
    TEST_CHECK(Rrms(db2)==approx(3.78810001,1e-6));

    DensityType dbn;
    TEST_CHECK( !(!db1));
    TEST_CHECK( !dbn);
    TEST_CHECK( !(!db1));
    TEST_CHECK( !dbn);
    }
