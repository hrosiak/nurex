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

TEST_CASE("smooth"){
    for(int a = 58; a<70; a++){
      double r1 = fusion_barrier(2,1,a,57);
      double r2 = fusion_barrier(2,1,a+1,57);
      CHECK(std::abs(r1-r2) < 0.5);
      CHECK(std::abs(fusion_barrier(1,1,a,57)-fusion_barrier(1,1,a+1,57)) < 0.5);
    }
}

TEST_CASE("barrier"){
  int A1, A2,Z1,Z2;
  double r;
  auto fx = [&](double x)->double{
            return nuclear_potential(x,A1,Z1,A2,Z2);
        };

  A1 = 12;
  Z1 = 6;
  A2 = 4;
  Z2 = 2;
  r = golden_search(fx,0,20);
  CHECK( r == approx(7.1,1.0));
  CHECK(fusion_barrier(A1,Z1,A2,Z2) == approx(2.35,0.5));

  A1 = 40;
  Z1 = 20;
  A2 = 4;
  Z2 = 2;
  r = golden_search(fx,0,20);
  CHECK( r == approx(7.7,1.0));
  CHECK(fusion_barrier(A1,Z1,A2,Z2) == approx(6.8,0.75));
  CHECK(fusion_barrier(44,20,4,2) == approx(6.4,0.75));
  CHECK(fusion_barrier(44,20,4,2) < fusion_barrier(40,20,4,2));

  CHECK(fusion_barrier(48,22,4,2) == approx(7.1,0.7));
  CHECK(fusion_barrier(51,23,4,2) == approx(7.4,0.7));

  CHECK(fusion_barrier(59,27,4,2) == approx(8.3,0.7));

  A1 = 162;
  Z1 = 66;
  A2 = 4;
  Z2 = 2;
  r = golden_search(fx,0,20);
  CHECK( r == approx(10.3,1.0));
  CHECK(fusion_barrier(A1,Z1,A2,Z2) == approx(17.1,1));

  CHECK(fusion_barrier(208,82,4,2) == approx(20.3,0.7));

  A1 = 209;
  Z1 = 89;
  A2 = 4;
  Z2 = 2;
  r = golden_search(fx,0,20);
  CHECK( r == approx(10.0,1.0));
  CHECK(fusion_barrier(A1,Z1,A2,Z2) == approx(21.0,1));

  CHECK(fusion_barrier(233,92,4,2) == approx(22.2,0.7));
}

TEST_CASE("heavy"){
    CHECK(fusion_barrier(16,8,58,28) == approx(31.8,1.0));
    CHECK(fusion_barrier(16,8,60,28) == approx(31.6,1.0));

    CHECK(fusion_barrier(40,20,124,50) == approx(118.7,3.0));
    CHECK(fusion_barrier(40,20,132,50) == approx(117.2,3.0));
}

TEST_CASE("proton"){
    CHECK(fusion_barrier(51,23,1,1) == approx(4.2,0.7));
    CHECK(fusion_barrier(59,27,1,1) == approx(5.0,0.7));
    CHECK(fusion_barrier(63,29,1,1) == approx(5.0,0.7));
    CHECK(fusion_barrier(65,29,1,1) == approx(4.9,0.7));
    CHECK(fusion_barrier(92,40,1,1) == approx(6.4,0.7));
    CHECK(fusion_barrier(94,40,1,1) == approx(6.2,0.7));
    CHECK(fusion_barrier(93,41,1,1) == approx(6.5,0.7));
    CHECK(fusion_barrier(95,42,1,1) == approx(6.65,0.7));
    CHECK(fusion_barrier(98,42,1,1) == approx(6.4,0.7));
    CHECK(fusion_barrier(103,45,1,1) == approx(6.7,0.7));
    CHECK(fusion_barrier(110,46,1,1) == approx(6.8,0.7));
    CHECK(fusion_barrier(109,47,1,1) == approx(7.0,0.7));
    CHECK(fusion_barrier(115,50,1,1) == approx(7.5,0.7));
    CHECK(fusion_barrier(124,50,1,1) == approx(7.3,0.7));
    
    CHECK(fusion_barrier(94,40,1,1) < fusion_barrier(92,40,1,1));
    CHECK(fusion_barrier(98,42,1,1) < fusion_barrier(95,42,1,1));

}

TEST_CASE("parametrized proton"){
    CHECK(fusion_barrier_parametrized(51,23,1,1) == approx(4.2,0.7));
    CHECK(fusion_barrier_parametrized(59,27,1,1) == approx(4.8,0.7));
    CHECK(fusion_barrier_parametrized(63,29,1,1) == approx(5.0,0.7));
    CHECK(fusion_barrier_parametrized(65,29,1,1) == approx(4.9,0.7));
    CHECK(fusion_barrier_parametrized(92,40,1,1) == approx(6.4,0.7));
    CHECK(fusion_barrier_parametrized(94,40,1,1) == approx(6.2,0.7));
    CHECK(fusion_barrier_parametrized(93,41,1,1) == approx(6.5,0.7));
    CHECK(fusion_barrier_parametrized(95,42,1,1) == approx(6.65,0.7));
    CHECK(fusion_barrier_parametrized(98,42,1,1) == approx(6.4,0.7));
    CHECK(fusion_barrier_parametrized(103,45,1,1) == approx(6.7,0.7));
    CHECK(fusion_barrier_parametrized(110,46,1,1) == approx(6.8,0.7));
    CHECK(fusion_barrier_parametrized(109,47,1,1) == approx(7.0,0.7));
    CHECK(fusion_barrier_parametrized(115,50,1,1) == approx(7.5,0.7));
    CHECK(fusion_barrier_parametrized(124,50,1,1) == approx(7.3,0.7));
    
    CHECK(fusion_barrier_parametrized(94,40,1,1) < fusion_barrier_parametrized(92,40,1,1));
    CHECK(fusion_barrier_parametrized(98,42,1,1) < fusion_barrier_parametrized(95,42,1,1));

}

TEST_CASE("parametrized alpha"){
  CHECK(fusion_barrier_parametrized(12,6,4,2) == approx(2.35,0.5));
  CHECK(fusion_barrier_parametrized(40,20,4,2) == approx(6.8,0.75));
  CHECK(fusion_barrier_parametrized(44,20,4,2) == approx(6.4,0.75));
  CHECK(fusion_barrier_parametrized(44,20,4,2) < fusion_barrier_parametrized(40,20,4,2));
  CHECK(fusion_barrier_parametrized(48,22,4,2) == approx(7.1,0.7));
  CHECK(fusion_barrier_parametrized(51,23,4,2) == approx(7.4,0.7));
  CHECK(fusion_barrier_parametrized(59,27,4,2) == approx(8.3,0.7));
  CHECK(fusion_barrier_parametrized(208,82,4,2) == approx(20.3,1.0));
  CHECK(fusion_barrier_parametrized(233,92,4,2) == approx(22.2,1.0));
}

TEST_CASE("parametrized heavy"){
    CHECK(fusion_barrier_parametrized(16,8,58,28) == approx(31.8,1.0));
    CHECK(fusion_barrier_parametrized(16,8,60,28) == approx(31.6,1.0));

    CHECK(fusion_barrier_parametrized(40,20,124,50) == approx(118.7,3.0));
    CHECK(fusion_barrier_parametrized(40,20,132,50) == approx(117.2,3.0));
}


TEST_CASE("par vs calculated"){
    std::vector<std::pair<int,int>> targets = {{16,8},{40,20},{124,50},{238,92}};
    for(auto     e : targets){
            CHECK(fusion_barrier_parametrized(e.first,e.second,2,1) == approx(fusion_barrier(e.first,e.second,2,1),1));
            CHECK(fusion_barrier_parametrized(e.first,e.second,3,1) == approx(fusion_barrier(e.first,e.second,3,1),1));
            CHECK(fusion_barrier_parametrized(e.first,e.second,3,2) == approx(fusion_barrier(e.first,e.second,3,2),1));
    }
}
