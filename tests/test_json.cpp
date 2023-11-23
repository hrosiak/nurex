#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include <vector>
#include <fstream>
#include "nurex/json_utils.h"

using namespace nurex;

json load_json(const char *fname){
    std::vector<std::string> res;
    std::ifstream jfile(fname,std::ifstream::in);

    if(!jfile){
        throw std::invalid_argument("Could not open config file");
    }

    std::string content;
    jfile.seekg(0, std::ios::end);
    content.resize(jfile.tellg());
    jfile.seekg(0, std::ios::beg);
    jfile.read(&content[0], content.size());
    jfile.close();

    json j = string_json(content);
    return j;
        //cout<<"JSON parsing error\n";
        //throw std::invalid_argument("Could not parse json file");

};
const std::string c_wrong=R"({
               "model":"XXX",
               "target":"12c",
               "projectile": "12c,
               })";

const std::string c1=R"({
               "model":"OLAFR",
               "target":{
                        "nucleus":"12c",
                        "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                        "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                       },

               "projectile":{
                           "nucleus":"64ni",
                           "proton_density":{"type":"fermi","parameters":[4.2413494,0.50732378]},
                           "neutron_density":{"type":"fermi","parameters":[4.2413494,0.50732378]}
                          }
               })";

const std::string c2=R"({
               "model":"OLAFR_FM",
               "target":{
                        "nucleus":"12c",
                        "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                        "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                       },

               "projectile":{
                           "nucleus":"64ni",
                           "proton_density":{"type":"fermi","parameters":[4.2413494,0.50732378]},
                           "neutron_density":{"type":"fermi","parameters":[4.2413494,0.50732378]}
                          }
               })";


const std::string c_12c_12c=R"({
               "model":"OLAFR_FM",
                "charge_changing_correction":"PRC82",
               "target":{
                        "nucleus":"12c",
                        "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                        "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                       },

               "projectile":{
                            "nucleus":"12c",
                            "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                            "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                          }
               })";

const std::string c_12c_12c_cc2=R"({
               "model":"OLAFR_FM",
                "charge_changing_correction":"PRC82",
                "coulomb_correction":"sommerfeld",
               "target":{
                        "nucleus":"12c",
                        "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                        "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                       },

               "projectile":{
                            "nucleus":"12c",
                            "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                            "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                          }
               })";

const std::string cexc=R"({
               "model":"OLAFR_FM",
                "charge_changing_correction":"evaporation",
                "coulomb_correction":"sommerfeld",
               "target":{
                        "nucleus":"12c",
                        "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                        "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                       },

               "projectile":{
                            "nucleus":"12c",
                            "proton_density":{"type":"ho","parameters":[1.548415436,1.6038565]},
                            "neutron_density":{"type":"ho","parameters":[1.548415436,1.6038565]}
                          },
                "excitation_function":20,
                "pairing_backshift":1
               })";

    TEST_CASE("json model"){
        json j;
        try{
            j = load_json("../../tests/config1.js");
        }
        catch(...){
            j = load_json("../tests/config1.js");
        }
        auto gm = json_model(j);
        CHECK_FALSE(!gm);
        CHECK_FALSE(!gm.Projectile());
        CHECK_FALSE(!gm.Target());
        CHECK(gm.Target().A()==12);
        CHECK(gm.Target().Z()==6);

        CHECK(gm.SigmaR(100)==approx(2050,10));
        CHECK(gm.SigmaR(200)==approx(1806,10));
    }
    
    TEST_CASE("json model from string"){
        auto gm = json_model(c1);
        CHECK_FALSE(!gm);
        CHECK(gm.Target().A()==12);
        CHECK(gm.Target().Z()==6);
        CHECK(gm.Projectile().A()==64);
        CHECK(gm.Projectile().Z()==28);
        CHECK(gm.SigmaR(30)==approx(2392,10));
        CHECK(gm.SigmaR(100)==approx(1943,10));
        CHECK(gm.SigmaR(200)==approx(1805,10));

        auto gm2 = json_model(c2);
        CHECK_FALSE(!gm2);
        CHECK(gm2.Target().A()==12);
        CHECK(gm2.Target().Z()==6);
        CHECK(gm2.Projectile().A()==64);
        CHECK(gm2.Projectile().Z()==28);
        CHECK(gm2.SigmaR(30)==approx(2745,10));
        CHECK(gm2.SigmaR(100)==approx(2106,10));
        CHECK(gm2.SigmaR(200)==approx(1835,10));

        gm2.SetCoulombCorrection(coulomb_correction_t::classic);
        CHECK(gm2.SigmaR(30)==approx(0.9139*2745,10));
        CHECK(gm2.SigmaR(100)==approx(0.97097*2106,10));
    }
    
    TEST_CASE("cc correction from json"){
        auto gmc = json_model(c_12c_12c);
        CHECK_FALSE(!gmc);
        CHECK(SigmaCC(gmc,30)==approx(1265.9*nurex::sigma_cc_scaling_factor(30),6));
        CHECK(SigmaCC(gmc,100)==approx(867.0*nurex::sigma_cc_scaling_factor(100),4));
        CHECK(SigmaCC(gmc,200)==approx(689.7*nurex::sigma_cc_scaling_factor(200),3));

    }
    
    TEST_CASE("couloumb correction"){
        auto gm1 = json_model(c_12c_12c);
        auto gm2 = json_model(c_12c_12c_cc2);
        CHECK(SigmaCC(gm1,100)>SigmaCC(gm2,100));
        CHECK(SigmaCC(gm1,50)>SigmaCC(gm2,50));
        CHECK(SigmaCC(gm1,2000)==approx(SigmaCC(gm2,2000),10));
    }
    
    TEST_CASE("wrong json config"){
        auto gm = json_model(c_wrong);
        CHECK(!gm);
    }
    
    TEST_CASE("DenityTable from json"){
        DensityFermi df(4,0.5);
        double r = Rrms(df);
        const int n = 100;
        auto b = linspace_array<n>(0,5*r);
        auto b2 = linspace_array<40>(0,4*r);
        std::vector<double> x,y;
        x.reserve(n);y.reserve(n);
        json j;
        j["type"] = "table";
        json pars;
        for(auto _b:b){
            json data;
            data.push_back(_b);
            data.push_back(df.Density(_b));
            pars.push_back(data);
        }
        j["parameters"] = pars;
        auto dt = json_density(j);
        CHECK_FALSE(!dt);
        for(auto _b:b2){
            CHECK(dt.Density(_b)==approx(df.Density(_b)).R(1e-3));
        }
        CHECK(Rrms(dt) == approx(Rrms(df),1e-4));
    }
    
    TEST_CASE("excitation parameters"){
        auto gm1 = json_model(cexc);
        CHECK_FALSE(!gm1);
    }
