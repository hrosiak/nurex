#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include <vector>
#include "nurex/NucID.h"

using namespace nurex;


    TEST_CASE("nucid"){
       CHECK(nucleus_id(1,1)==10010);
       CHECK(nucleus_id(5,1)==50010);
       CHECK(nucleus_id(12,6)==120060);
    }
    
    TEST_CASE("nucid inverse"){
       CHECK(nucleus_from_id(10010)[0]==1);
       CHECK(nucleus_from_id(10010)[1]==1);
       CHECK(nucleus_from_id(50010)[0]==5);
       CHECK(nucleus_from_id(10010)[1]==1);
       CHECK(nucleus_from_id(120060)[0]==12);
       CHECK(nucleus_from_id(120060)[1]==6);
    }
    
    TEST_CASE("nuc from symbol "){
       auto r = nucleus_from_symbol("12C");
       CHECK(r[0]==12);
       CHECK(r[1]==6);
       r = nucleus_from_symbol("12c");
       CHECK(r[0]==12);
       CHECK(r[1]==6);

       r = nucleus_from_symbol("238U");
       CHECK(r[0]==238);
       CHECK(r[1]==92);

       r = nucleus_from_symbol("122Rx");
       CHECK(r[0]==0);
       CHECK(r[1]==0);       
    }
    
    TEST_CASE("symbol from nuc"){       
       CHECK(nucleus_symbol(12,6).compare("12C")==0);    
       CHECK(nucleus_symbol(0,6).compare("")==0);
       CHECK(nucleus_symbol(48,20).compare("48Ca")==0);    
    }
