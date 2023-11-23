#include <math.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include "nurex/nurex.h"
#include "nurex/json_utils.h"
#include "json.hpp"

using namespace std;
using namespace nurex;
using json = nlohmann::json;

void help(){
        std::cout<<"usage: nurex_calculator [options]";
        std::cout<<"options:\n";
        std::cout<<"\t--model MODEL_NAME - ie OLAZR, MOLFR_FM ...\n";
        std::cout<<"\t--range value - optional to change range parameter, default 0.39\n";
        std::cout<<"\t--target - target symbol, ie 12c\n";
        std::cout<<"\t--projectile - projectile symbol, ie 27al\n";
        std::cout<<"\t--energy - energy in Mev/u, for range type Emin,Emax,Estep\n";
        std::cout<<"\t--file FILENAME - save results to file FILENAME\n";
        std::cout<<"\t--config FILENAME - load config from JSON formatted file";
        std::cout<<"example:\n";
        std::cout<<"nurex_calculator --model MOLZR_FM --target 12c --projectile 12c --energy 400\n";
        std::cout<<"nurex_calculator --model OLAFR --target 12c --projectile 12c --energy 100,1000,10\n";
        std::cout<<"nurex_calculator --config config.json\n";
}

void print_density(const DensityType&);
json load_json(const char *fname);
char* getCmdOption(char ** begin, char ** end, const std::string & option);

int main( int argc, char * argv[] )
{
    std::string energy_argument;
    std::string filename;
    std::vector<double> energies;
    bool fileout = false;
    Nucleus projectile;
    Nucleus target;
    double custom_range = 0.0;
    json j;

    if(argc == 1 ){
        help();
        return 0;
    }

    auto arg = getCmdOption(argv,argv+argc,"--config");
    if(argc==2 || arg){
        try{
            char *config_file = (argc==2)?argv[1]:arg;
            j = load_json(config_file);
            energies = json_energies(j);
            if(j.count("file")){
                    if(j["file"].is_string())filename = j["file"].get<std::string>();
                    fileout = true;
                }
        }
        catch(...){
            cout<<"Could not load the json config file"<<"\n";
            return 0;
        }
    }
    else{ // settings from command line arguments
        arg = getCmdOption(argv,argv+argc,"--model");
        if(arg){
            j["model"] = arg;
            }

        arg = getCmdOption(argv,argv+argc,"--range");
        if(arg){
            std::string val(arg);
            j["range"] = std::stod(val);
            }

        arg = getCmdOption(argv,argv+argc,"--target");
        if(arg){
            j["target"] = arg;
            }

        arg = getCmdOption(argv,argv+argc,"--projectile");
        if(arg){
            j["projectile"] = arg;
            }

        arg = getCmdOption(argv,argv+argc,"--file");
        if(arg){
            j["file"] = arg;
        }

        arg = getCmdOption(argv,argv+argc,"--energy");
        if(arg){
            energy_argument = arg;
            double _e = atof(energy_argument.c_str());
            if(_e<8){
                printf("wrong energy value\n");
                return 0;
                }

            if(energy_argument.find(",")==std::string::npos){
                energies.push_back(_e);
            }
            else{
                double emin = 0, emax = 1500, estep = 50;
                auto pos1 = energy_argument.find(",");
                auto s1 = energy_argument.substr(0,pos1);
                auto s2 = energy_argument.substr(pos1+1);

                _e = atof(s1.c_str());
                emin=_e;
                _e = atof(s2.c_str());
                emax=_e;
                estep=atof(s2.substr(s2.find(",")+1).c_str());

                if(emin<8 || emax<8 || emax<=emin || estep<=0){
                    printf("wrong energy value\n");
                    return 0;
                }

            for(double e = emin;e<=emax;e+=estep){
                energies.push_back(e);
                }
            }
            }

        }//end of else json config

    GlauberModelType gm;
    try{
        gm = json_model(j);
    }
    catch(...){
        cout<<"Could not parse json config file. Wrong parameters or json format."<<"\n";
        return 0;
        }

    if(!gm){
        printf("Wrong json config file or arguments\n");
        return 0;
    }

    std::cout<<"Model = "<<j["model"]<<"\n";
    if(custom_range>0.0){
        std::cout<<"range parameter: "<<custom_range<<" fm\n";
    }
    std::cout<<"--- Projectile: "<<symbol(gm.Projectile())<<" ---\n";
    std::cout<<"protons: Rrms = "<<gm.Projectile().RrmsProton()<<"fm, ";
    print_density(gm.Projectile().GetDensityProton());
    std::cout<<"\n";
    std::cout<<"neutrons: Rrms = "<<gm.Projectile().RrmsNeutron()<<"fm, ";
    print_density(gm.Projectile().GetDensityNeutron());

    std::cout<<"\n--- Target: "<<symbol(gm.Target())<<" ---\n";
    std::cout<<"protons: Rrms = "<<gm.Target().RrmsProton()<<"fm, ";
    print_density(gm.Target().GetDensityProton());
    std::cout<<"\n";
    std::cout<<"neutrons: Rrms = "<<gm.Target().RrmsNeutron()<<"fm, ";
    print_density(gm.Target().GetDensityNeutron());

    ofstream fw;
    if(fileout){
        fw.open(filename.c_str());
        if(!fw.is_open()){
            printf("Could not open file\n");
            return 0;
        };
    }

    double sigma_r, sigma_cc;

    std::cout<<"\n--- Cross Sections: ---\n";
    if(energies.size()==0){
        cout<<"no energy specified\n";
        return 0;
    }
    for(double e:energies){
        sigma_r = gm.SigmaR(e);
        sigma_cc= gm.SigmaCC(e);
        if(fileout){
            fw<<e<<" "<<sigma_r<<" "<<sigma_cc<<" ";
          //  fw<<gm->CoulombCorrection(gm->SigmaR_prev)<<" ";
//            fw<<sigma_cc_scaling_factor(e);
            fw<<"\n";
        }
        std::cout<<"Energy = "<<e<<" MeV/u, Sigma_R = "<<sigma_r<<" mb, Sigma_CC = "<<sigma_cc<<" mb";
        //std::cout<<", Coulomb Correction = "<<gm->CoulombCorrection(gm->SigmaR_prev);
        //std::cout<<", CC energy scaling = "<<sigma_cc_scaling_factor(e);
        std::cout<<"\n";
    }

    if(fileout){
        fw.close();
    }

    return 1;
}

void print_density(const DensityType& df){
    if(df.type()==density_type::fermi){
        cout<<"density type: Fermi";
        cout<<", parameters: "<<df.GetParameter(0)<<", "<<df.GetParameter(1);
        if(df.GetParameter(2)>0.0){
          cout<<", "<<df.GetParameter(2);
        }
    }

    if(df.type()==density_type::ho){
        cout<<"density type: HO";
        cout<<", parameters: "<<df.GetParameter(0)<<", "<<df.GetParameter(1);
    }

    if(df.type()==density_type::zero){
        cout<<"density type: Zero";
    }

    if(df.type()==density_type::dirac){
        cout<<"density type: Dirac";
    }

    if(df.type()==density_type::table){
        cout<<"density type: Table";
    }
}

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

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return nullptr;
}
