#include <nlohmann/json.hpp>
#include "nurex/json_utils.h"
#include "nurex/numerics.h"
#include "nurex/Models.h"
#include "nurex/DefaultDensities.h"
using json = nlohmann::json;

namespace nurex{

bool json_test(const json &j){
    if( !j.count("model")|| !j.count("projectile") || !j.count("target")){
        return false;
    }
    else{
        return true;
    }
};

json string_json(const std::string& string){
    try{
        auto j = json::parse(string);
        return j;
    }
    catch(...){
        json j;
        return j;
    }
    };


DensityType json_density(const json &j){
    std::string density;
    double norm = 1.0;
    if(j.count("norm")){
	norm = j["norm"].get<double>();
    }
    if(j.count("type") && j.count("parameters") && j["type"].is_string()){
        density = j["type"].get<std::string>();

        if( density.compare("Fermi")==0
          || density.compare("fermi")==0){
            if(j["parameters"].is_array()){
                double p0 = j["parameters"].at(0).get<double>();
                double p1 = j["parameters"].at(1).get<double>();
                double p2 = 0.0;
                if(j["parameters"].size()>2)p2 = j["parameters"].at(2).get<double>();
                return DensityType(make_density<DensityFermi>(p0, p1, p2, norm));
            }

        }

        if( density.compare("HO")==0 || density.compare("ho")==0 ){
            if(j["parameters"].is_array())
                return DensityType(make_density<DensityHO>(j["parameters"].at(0).get<double>(),j["parameters"].at(1).get<double>(), norm));
        }
        if( density.compare("Dirac")==0 || density.compare("dirac")==0 ){
                return DensityType(make_density<DensityDirac>(norm));
        }
        if( density.compare("zero")==0 || density.compare("Zero")==0 ){
                return DensityType(make_density<DensityZero>());
        }
        if( density.compare("gaussian")==0 || density.compare("Gaussian")==0 ){
                return DensityType(make_density<DensityGaussian>(j["parameters"].at(0).get<double>(), norm));
        }
        if( density.compare("table")==0 || density.compare("Table")==0 ){
                if(j["parameters"].is_array() && j["parameters"].size()>1){
                    std::vector<double> b;
                    std::vector<double> d;
                    const auto num = j["parameters"].size();
                    b.reserve(num);
                    d.reserve(num);
                    for(auto& e:j["parameters"]){
                        if(e.size()!=2 || !e[0].is_number() || !e[0].is_number())return DensityType();
                        b.push_back(e.at(0).get<double>());
                        d.push_back(e.at(1).get<double>());
                    }
                    return DensityType(make_density<DensityTable>(b,d, norm));
                }
        }

#ifndef NO_FILESYSTEM
        if( density.compare("file")==0 || density.compare("File")==0 ){
            if(j["parameters"].is_string()){
                return DensityType(density_from_file(j["parameters"].get<std::string>().c_str()));
            }
        }
#endif
    }
    return DensityType();
}

std::vector<double> json_energies(const json &j){
    std::vector<double> energies;
    if(j.count("energy")){
                auto e = j.at("energy");
                if(e.is_number()){
                    energies.push_back(j["energy"].get<double>());
                    }
                if(e.is_string()){
                    double _e = std::stod(j["energy"].get<std::string>());
                    energies.push_back(_e);
                    }
                if(e.is_array()){
                    for(auto &el:e){
                        if(el.is_number())
                        energies.push_back(el.get<double>());
                        }
                    }
                if(e.is_object()){
                    if(e.count("min")>0 && e.count("max")>0 && (e.count("num")>0 || e.count("step")>0)){
                        double emin = e["min"].get<double>();
                        double emax = e["max"].get<double>();
                        unsigned int num=0;
                        if(e.count("step")){
                            num = 1.0+(emax-emin)/e["step"].get<unsigned int>();
                        }
                        if(e.count("num")){
                            num = e["num"].get<unsigned int>();
                        }
                        energies  =  nurex::linspace_vector(emin,emax,num);
                    }
                    }
    }//end of json energy check
    return energies;
}


Nucleus json_nucleus(const json& j){
    std::string symbol;
    DensityType dp,dn;

    if(j.is_string()){
        symbol = j.get<std::string>();
        return  get_default_nucleus(symbol.c_str());
    }

    if(j.is_array() && j.size()==2){
        return  get_default_nucleus(j[0].get<int>(), j[1].get<int>());
    }

    if(j.count("nucleus")){
        if(j["nucleus"].is_string())symbol = j["nucleus"].get<std::string>();
        }
    if(j.count("proton_density")){
        auto obj = j["proton_density"];
        dp = json_density(obj);
    }
    if(j.count("neutron_density")){
        auto obj = j["neutron_density"];
        dn = json_density(obj);
    }

    if(!dp || !dn){
        return  Nucleus();
    }
    else{
        auto nuc = nucleus_from_symbol(symbol);
        return Nucleus(nuc[0],nuc[1],dp,dn);
    }
}

EvaporationParameters json_evaporation_parameters(const json &j){
    EvaporationParameters evaporation_parameters;
    if(j.count("excitation_energy")){
        evaporation_parameters.Emax = j["excitation_energy"].get<double>();
    }
    // neutron removal cross sections
    if(j.count("xn_cs")){
           if(j["xn_cs"].is_string()){
            std::string type = j["xn_cs"].get<std::string>();
            if(type.compare("epax") == 0){evaporation_parameters.xn_cs = 1;}
        }
        else {
           evaporation_parameters.xn_cs = j["xn_cs"].get<int>();
        }
    }

    if(j.count("evaporation_preset") && j["evaporation_preset"].is_string()){
        std::string type = j["evaporation_preset"].get<std::string>();
        if(type.compare("abla") == 0){
            evaporation_parameters.preset = evaporation_preset_type::abla;
            }
        if(type.compare("nurex") == 0){
            evaporation_parameters.preset = evaporation_preset_type::nurex;
            }
    }

    if(j.count("level_density") && j["level_density"].is_string()){
        std::string type = j["level_density"].get<std::string>();
        if(type.compare("GC_GEM") == 0){
            evaporation_parameters.density = level_density_type::GC_GEM;
            }
        else if(type.compare("GC_RIPL") == 0){
            evaporation_parameters.density = level_density_type::GC_RIPL;
            }
        else if(type.compare("GC_KTUY05") == 0){
            evaporation_parameters.density = level_density_type::GC_KTUY05;
            }
        else if(type.compare("ABLA") == 0){
            evaporation_parameters.density = level_density_type::ABLA;
            }
    }

    if(j.count("coulomb_barrier") && j["coulomb_barrier"].is_string()){
        std::string type = j["coulomb_barrier"].get<std::string>();

        if(type.compare("bass") == 0 || type.compare("bass80") == 0){
            evaporation_parameters.barrier = barrier_type::bass80;
            }
        if(type.compare("none") == 0){
            evaporation_parameters.barrier = barrier_type::none;
            }
        else if (type.compare("parametrized") == 0) {
            evaporation_parameters.barrier = barrier_type::parametrized;
            }
        else evaporation_parameters.barrier = 0;
    }


    // eveporation configuration
    if(j.count("disable_neutron_evaporation") && j["disable_neutron_evaporation"].get<bool>()){
        evaporation_parameters.config = evaporation_parameters.config
                                       |evaporation_config_type::disable_neutron;
    }
    if(j.count("disable_imf_evaporation") && j["disable_imf_evaporation"].get<bool>()){
        evaporation_parameters.config = evaporation_parameters.config
                                       |evaporation_config_type::disable_imf;
    }

    return evaporation_parameters;
}

GlauberModelType json_model(const json &j){
    std::string model;
    GlauberModelType gm;
    using fmtype = NNCrossSection_FermiMotion;
    range_parameters custom_range;
    Nucleus target;
    Nucleus projectile;
    double fmp1=0.0,fmp2=0.0;
    double fc=0.0;
    EvaporationParameters evaporation_parameters;
    
    if(j.count("model"))model = j["model"];
    if(j.count("range") && j["range"].is_number_float()){
        custom_range.pp = j["range"].get<double>();
        custom_range.pn = custom_range.pp;
        }
    else if(j.count("range") && j["range"].is_array() && (j["range"].size()>1)){
        custom_range.pp = j["range"][0].get<double>();
        custom_range.pn = j["range"][1].get<double>();
    }
    if(j.count("fermi_motion") && j["fermi_motion"].is_array()){
        fmp1 = j["fermi_motion"][0].get<double>();
        if(fmp1<0)fmp1=0;
        fmp2 = j["fermi_motion"][1].get<double>();
        if(fmp2<0)fmp2=0;
        }
    if(j.count("fermi_energy_coefficient")){
        fc = j["fermi_energy_coefficient"].get<double>();
        if(fc<0)fc=0.0;
        }
    
    if(j.count("target")){
        target = json_nucleus(j["target"]);
        }
    
    if(j.count("projectile")){
        projectile = json_nucleus(j["projectile"]);
        }

    if(!target || !projectile){        
        assert(((void)"target or projectile is wrong", false));
        return gm;
        }

    if( (model.compare("MOL")==0) || (model.compare("MOLFR")==0) || (model.compare("MOLZR")==0)){
        using MT = GlauberModel<MOL>;
        gm = make_glauber_model<MOL>(projectile,target);
    }
    else if( (model.compare("MOLZR_FM")==0) || (model.compare("MOL_FM")==0) || (model.compare("MOLFR_FM")==0)){
        using MT = GlauberModel<MOL, fmtype>;
        gm = make_glauber_model<MOL,fmtype>(projectile,target);
        if(projectile.A()==1 ||target.A()==1){
                static_cast<MT*>(gm.get_object())->get_sigma_nn().SetMomentum(0, 90.0);
                }
        if(fmp1>0.0 || fmp2>0.0)static_cast<MT*>(gm.get_object())->get_sigma_nn().SetMomentum(fmp1, fmp2);
    }
    else if( (model.compare("MOL4C")==0) || (model.compare("MOL4CFR")==0) || (model.compare("MOL4CZR")==0)){
        using MT = GlauberModel<MOL4C>;
        gm = make_glauber_model<MOL4C>(projectile,target);
    }
    else if( (model.compare("MOL4CZR_FM")==0) || (model.compare("MOL4C_FM")==0) || (model.compare("MOL4CFR_FM")==0)){
        using MT = GlauberModel<MOL4C, fmtype>;
        gm = make_glauber_model<MOL4C,fmtype>(projectile,target);
        if(projectile.A()==1 ||target.A()==1){
                static_cast<MT*>(gm.get_object())->get_sigma_nn().SetMomentum(0, 90.0);
                }
        if(fmp1>0.0 || fmp2>0.0)static_cast<MT*>(gm.get_object())->get_sigma_nn().SetMomentum(fmp1, fmp2);
    }
    else if( (model.compare("OLA")==0) || (model.compare("OLAFR")==0) || (model.compare("OLAZR")==0)){
        using MT = GlauberModel<OLA>;
        gm = make_glauber_model<OLA>(projectile,target);
    }
    else if( (model.compare("OLAZR_FM")==0) || (model.compare("OLA_FM")==0) || (model.compare("OLAFR_FM")==0)){
        using MT = GlauberModel<OLA, fmtype>;
        gm = make_glauber_model<OLA,fmtype>(projectile,target);
        if(projectile.A()==1 ||target.A()==1){
                static_cast<MT*>(gm.get_object())->get_sigma_nn().SetMomentum(0, 90.0);
                }
        if(fmp1>0.0 || fmp2>0.0)static_cast<MT*>(gm.get_object())->get_sigma_nn().SetMomentum(fmp1, fmp2);
    }

    ////////// testing part
    else if( (model.compare("OLAZR_FMD")==0) || (model.compare("OLA_FMD")==0)|| (model.compare("OLAFR_FMD")==0)){
            using MT = GlauberModel<OLA_FMD,NNCrossSection_FermiMotionD>;
            gm = MT(projectile,target);
            if(fc>0.0){
                static_cast<MT*>(gm.get_object())->SetFECoefficient(fc);
            }
    }
    else if( (model.compare("MOL4CZR_FMD")==0) || (model.compare("MOL4CFR_FMD")==0)|| (model.compare("MOL4C_FMD")==0)){
            using MT = GlauberModel<MOL4C_FMD,NNCrossSection_FermiMotionD>;
            gm = MT(projectile,target);
            if(fc>0.0){
                static_cast<MT*>(gm.get_object())->SetFECoefficient(fc);
            }
    }
    else if( (model.compare("MOLZR_FMD")==0) || (model.compare("MOLFR_FMD")==0)|| (model.compare("MOL_FMD")==0)){
            using MT = GlauberModel<MOL_FMD,NNCrossSection_FermiMotionD>;
            gm = MT(projectile,target);
            if(fc>0.0){
                static_cast<MT*>(gm.get_object())->SetFECoefficient(fc);
            }
    }
    else{
        assert(((void)"glauber model id unknown", false));
        return gm;
    }
    // finite range models so check  range parameters
    if(   (model.compare("OLAFR_FMD")==0)
        ||(model.compare("OLAFR_FM")==0)
        ||(model.compare("OLAFR")==0)                
        ||(model.compare("MOL4CFR_FMD")==0)
        ||(model.compare("MOL4CFR_FM")==0)
        ||(model.compare("MOL4CFR")==0)
        ||(model.compare("MOLFR_FM")==0)
        ||(model.compare("MOLFR_FMD")==0)
        ){
        if(custom_range.is_zero()){  // if zero, set default range
            custom_range.pp = 0.39;
            custom_range.pn = 0.39;
            }
        }
    if(!custom_range.is_zero()){
            gm.SetRange(custom_range.pp, custom_range.pn);
        }

    /// seting coulomb correction
    if(j.count("coulomb_correction") && j["coulomb_correction"].is_string()){
        coulomb_correction_t cc_type = coulomb_correction_t::none;
        std::string type = j["coulomb_correction"];
        if(type.compare("classic") == 0){cc_type = coulomb_correction_t::classic;}
        if(type.compare("relativistic") == 0){cc_type = coulomb_correction_t::relativistic;}
        if(type.compare("sommerfeld") == 0){cc_type = coulomb_correction_t::sommerfeld;}
        if(type.compare("none") == 0){cc_type = coulomb_correction_t::none;}
        gm.SetCoulombCorrection(cc_type);
    }

    /// seting charge-changing correction
    if(j.count("charge_changing_correction") && j["charge_changing_correction"].is_string()){
        cc_correction_t ccc_type = cc_correction_t::none;
        std::string type = j["charge_changing_correction"];
        if(type.compare("PRC82") == 0){ccc_type = cc_correction_t::PRC82;}
        if(type.compare("none") == 0){ccc_type = cc_correction_t::none;}
        if(type.compare("evaporation") == 0){ccc_type = cc_correction_t::evaporation;}
        gm.SetCCCorrection(ccc_type);
    }

    /// evaporation parameters
    if(j.count("excitation_energy")){
        evaporation_parameters.Emax = j["excitation_energy"].get<double>();
    }

    // neutron removal scaling
    if(j.count("n_removal_scaling")){
        evaporation_parameters.n_removal_scaling = j["n_removal_scaling"].get<double>();
    }

    // neutron removal cross sections
    if(j.count("xn_cs")){
           if(j["xn_cs"].is_string()){
            std::string type = j["xn_cs"].get<std::string>();
            if(type.compare("epax") == 0){evaporation_parameters.xn_cs = 1;}
        }
        else {
           evaporation_parameters.xn_cs = j["xn_cs"].get<int>();
        }
    }

    if(j.count("evaporation_preset") && j["evaporation_preset"].is_string()){
        std::string type = j["evaporation_preset"].get<std::string>();
        if(type.compare("abla") == 0){
            evaporation_parameters.preset = evaporation_preset_type::abla;
            }
        if(type.compare("nurex") == 0){
            evaporation_parameters.preset = evaporation_preset_type::nurex;
            }
    }

    if(j.count("level_density") && j["level_density"].is_string()){
        std::string type = j["level_density"].get<std::string>();
        if(type.compare("GC_GEM") == 0){
            evaporation_parameters.density = level_density_type::GC_GEM;
            }
        else if(type.compare("GC_RIPL") == 0){
            evaporation_parameters.density = level_density_type::GC_RIPL;
            }
        else if(type.compare("GC_KTUY05") == 0){
            evaporation_parameters.density = level_density_type::GC_KTUY05;
            }
        else if(type.compare("ABLA") == 0){
            evaporation_parameters.density = level_density_type::ABLA;
            }
    }

    if(j.count("coulomb_barrier") && j["coulomb_barrier"].is_string()){
        std::string type = j["coulomb_barrier"].get<std::string>();

        if(type.compare("bass") == 0 || type.compare("bass80") == 0){
            evaporation_parameters.barrier = barrier_type::bass80;
            }
        if(type.compare("none") == 0){
            evaporation_parameters.barrier = barrier_type::none;
            }
        else if (type.compare("parametrized") == 0) {
            evaporation_parameters.barrier = barrier_type::parametrized;
            }
        else evaporation_parameters.barrier = 0;
    }


    // eveporation configuration
    if(j.count("disable_neutron_evaporation")){
        if( (j["disable_neutron_evaporation"].is_boolean() && j["disable_neutron_evaporation"].get<bool>()) 
        || j["disable_neutron_evaporation"].is_number() && (j["disable_neutron_evaporation"].get<int>() > 0.0)
        ){
            evaporation_parameters.config = evaporation_parameters.config
                                       |evaporation_config_type::disable_neutron;
            } 
    }
    if(j.count("disable_imf_evaporation")){
        if(   (j["disable_imf_evaporation"].is_boolean() && j["disable_imf_evaporation"].get<bool>())
           || (j["disable_imf_evaporation"].is_number() && (j["disable_imf_evaporation"].get<int>()>0))
          ){
                evaporation_parameters.config = evaporation_parameters.config
                                       |evaporation_config_type::disable_imf;
            }        
    }

    //test
    if(j.count("evaporation_test")){
        evaporation_parameters.config = evaporation_parameters.config
                                       |evaporation_config_type::test;
    }

    gm.SetEvaporationParameters(evaporation_parameters);
    return gm;
}

GlauberModelType json_model(const std::string& string){
    return json_model(string_json(string));
};

std::vector<double> json_energies(const std::string& string){
    return json_energies(string_json(string));
};

DensityType json_density(const std::string& string){
    return json_density(string_json(string));
}

EvaporationParameters json_evaporation_parameters(const std::string& string){
    return json_evaporation_parameters(string_json(string));
}

}
