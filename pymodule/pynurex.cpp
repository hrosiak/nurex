#include <exception>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "nurex/nurex.h"
#include "nurex/evaporation.h"
#include "nurex/epax.h"
#include "nurex/ame16.h"
#include "nurex/json_utils.h"
//#include "nurex/data_evap.h"
#include <iostream>
namespace py = pybind11;
using namespace nurex;

GlauberModelType make_model(py::dict d){
    py::object dumps = py::module_::import("json").attr("dumps");
    std::string j = py::str(dumps(d));
    //std::string j = py::str(py::str(d).attr("replace")("'","\""));
    GlauberModelType gm{json_model(j)};
    if(!gm){
         throw std::invalid_argument("wrong model configuration");
    }
    return gm;
}

double py_sigma_r(GlauberModelType &m, double e){
    double res = m.SigmaR(e);
    assert(!isnan(res));
    return res;
}

double py_sigma_cc(GlauberModelType &m, double e){
    return m.SigmaCC(e);
}

DensityType make_density_from_dict(py::dict d){
    std::string j = py::str(py::str(d).attr("replace")("'","\""));
    return json_density(j);
}

py::object defaults_nucleus_python(py::args args){
    if(args.size()<2){
        Nucleus n = get_default_nucleus(std::string( py::str(args[0]) ).c_str());
        if(!n)return py::none();
        else return py::cast(n);
    }
    else{
        Nucleus n = get_default_nucleus( py::int_(args[0]), py::int_(args[1]));
        if(!n)return py::none();
        else return py::cast(n);
    }

}

EvaporationParameters make_evaporation_parameters(py::dict d){
    py::object dumps = py::module_::import("json").attr("dumps");
    std::string j = py::str(dumps(d));  
    return json_evaporation_parameters(j);  
}

/*
py::dict get_evaporation_data(int A, int Z){
    py::dict r;
    auto d = evaporation::get_data(A, Z);
    r["e"] = d.e;
    r["rn"] = d.rn;
    r["rch"] = d.rch;
    return r;
}
*/
py::dict density_object(DensityType const &d){
    py::dict r;
    py::list par;
    r["type"] = density_type_to_string(d);
    for(int i=0;i<d.NumberOfParameters();i++){
        par.append(d.GetParameter(i));
    }
    r["parameters"] = par;
    return r;
};

py::dict nucleus_object(Nucleus &n){
    py::dict r;
    py::list par;
    par.append(n.A());
    par.append(n.Z());
    r["nucleus"] = par;
    r["proton_density"] = density_object(n.GetDensityProton());
    r["neutron_density"] = density_object(n.GetDensityNeutron());
    return r;
};

double py_charge_evaporation_function(int A, int Z, double Ex, double j){
    prefragment n(A,Z);
    return charge_evaporation_function(n,Ex,j);
}

ExcitationFunction& get_custom_excitation(){
    return nurex::custom_excitation;
} 

void custom_excitation_set(std::vector<double>e, std::vector<double>w, int i){
    nurex::custom_excitation.set(e,w,i);
}

double custom_excitation_w(double e, int i){
    return nurex::custom_excitation.w(e,i);
}

void custom_excitation_reset(){
    nurex::custom_excitation.reset();
}

PYBIND11_MODULE(_pynurex,m){
     py::class_<DensityFermi>(m,"DensityFermi")
             .def(py::init<double, double, double, double>(), "constructor", py::arg("r"),py::arg("z"),py::arg("w")=0, py::arg("normalization")=1)
             .def("density",&DensityFermi::Density)
             .def("normalize",&DensityFermi::Normalize)
             .def("r_rms",[](const DensityFermi &d){return Rrms(d);});

     py::class_<DensityHO>(m,"DensityHO")
             .def(py::init<double, double, double>(), "constructor", py::arg("r"),py::arg("a"), py::arg("normalization")=1)
             .def("density",&DensityHO::Density)
             .def("normalize",&DensityHO::Normalize)
             .def("r_rms",[](const DensityHO &d){return Rrms(d);});

     py::class_<DensityDirac>(m,"DensityDirac")
             .def(py::init<double>(), "constructor", py::arg("normalization")=1)
             .def("density",&DensityDirac::Density)
             .def("normalize",&DensityDirac::Normalize)
             .def("r_rms",[](const DensityFermi &d){return 0.0;});

     py::class_<DensityZero>(m,"DensityZero")
             .def(py::init<>(), "constructor")
             .def("density",&DensityZero::Density)
             .def("normalize",&DensityZero::Normalize)
             .def("r_rms",[](const DensityFermi &d){return 0.0;});

     py::class_<DensityType>(m,"DensityType")
             //.def(py::init<DensityType&&>());
             .def(py::init<DensityFermi>())
             .def(py::init<DensityHO>())
             .def(py::init<DensityZero>())
             .def(py::init<DensityDirac>())
             .def("r_rms",&DensityType::Rrms,"returns Rrms of the density in fm")
             .def("density",&DensityType::Density,"density value at radius r", py::arg("r"))
             .def("number_of_parameters",&DensityType::NumberOfParameters)
             .def("get_parameter",&DensityType::GetParameter);

     py::class_<Nucleus>(m,"Nucleus")
             .def(py::init<int, int, DensityType, DensityType>(),"create nucleus by specifying A, Z and proton and neutron densities", py::arg("a"), py::arg("z"), py::arg("proton_density"), py::arg("neutron_density"))
             .def(py::init<const char*>(),"initialize nucleus specified by nucleus strin (ie \"12C\") using default density")
             .def("A",&Nucleus::A)
             .def("N",&Nucleus::N)
             .def("Z",&Nucleus::Z)
             .def("r_rms_proton",&Nucleus::RrmsProton,"returns Rrms of the proton density")
             .def("r_rms_neutron",&Nucleus::RrmsNeutron,"returns Rrms of the neutron density")
             .def("density_neutron",&Nucleus::DensityNeutron,"neutron density at radius",py::arg("r"))
             .def("density_proton", &Nucleus::DensityProton,"proton density at radius",py::arg("r"))
             .def("get_density_proton",&Nucleus::GetDensityProton)
             .def("get_density_neutron",&Nucleus::GetDensityNeutron);

     py::class_<NNCrossSectionFit>(m,"NNCrossSectionFit")
             .def(py::init<>(), "constructor")
             .def("np",&NNCrossSectionFit::np)
             .def("pp",&NNCrossSectionFit::pp);

     py::class_<GlauberModelType>(m,"GlauberModelType")
             .def("sigma_r",&GlauberModelType::SigmaR,"total reaction cross section", py::arg("energy"))
             .def("sigma_cc",&GlauberModelType::SigmaCC,"total charge changing cross section", py::arg("energy"))
             .def("sigma_xn",&GlauberModelType::SigmaXN,"total neutron removal cross section",py::arg("energy"))
             .def("sigma_ins",&GlauberModelType::SigmaINs, "neutron removals cross sectoin", py::arg("energy"))
             .def("set_evaporation_parameters",&GlauberModelType::SetEvaporationParameters)
             .def("get_evaporation_parameters",&GlauberModelType::GetEvaporationParameters)
             .def("set_excitation_function_type",&GlauberModelType::SetExcitationFunctionType)
             //.def("set_cc_correction",&GlauberModelType::SetCCCorrection)
             .def("set_range",&GlauberModelType::SetRange)
             .def("get_range",[](GlauberModelType& self){
                 py::dict r;
                 auto range = self.GetRange();
                 r["pp"] = range.pp;
                 r["pn"] = range.pn;
                 return r;
             })
             .def("X",&GlauberModelType::X)
             .def("Xpp",&GlauberModelType::Xpp)
             .def("Xpn",&GlauberModelType::Xpn)
             .def("Xnp",&GlauberModelType::Xnp)
             .def("Xnn",&GlauberModelType::Xnn)
             .def("projectile", &GlauberModelType::Projectile)
             .def("target", &GlauberModelType::Target)
	     .def("n_removals_evaporation",[](GlauberModelType& self){
			     py::dict r;
			     auto data = self.n_removals_evaporation( 0.0);
	             py::list ptot, pch, pn, pp, pd, pt, phe3, pa, pimf;
			     for(size_t i=0;i<data.Ptot.size();i++){
			     	ptot.append(data.Ptot[i]);
			    	pch.append(data.Pch[i]);
                    pn.append(data.Pn[i]);
                    pp.append(data.Pp[i]);
                    pd.append(data.Pd[i]);
                    pt.append(data.Pt[i]);
                    phe3.append(data.Phe3[i]);
                    pa.append(data.Pa[i]);
                    pimf.append(data.Pimf[i]);
			     	}
			     r["Ptot"] = ptot;
			     r["Pch"] = pch;
                 r["Pn"] = pn;
                 r["Pp"] = pp;
                 r["Pd"] = pd;
                 r["Pt"] = pt;
                 r["Phe3"] = phe3;
                 r["Pa"] = pa;
                 r["Pimf"] = pimf;
			     return r;
			     });

    py::class_<ExcitationFunction>(m,"ExcitationFunction")
            .def("set", &ExcitationFunction::set,"set excitation funtion for i hole",py::arg("energies"),py::arg("w"), py::arg("i"))
            .def("w", &ExcitationFunction::w,"get excitation funtion for excitation energy and i number of holes", py::arg("energy"), py::arg("i"))
            .def("emax", &ExcitationFunction::emax)
            .def("reset", &ExcitationFunction::reset);

    py::class_<prefragment>(m,"prefragment")
            .def(py::init<int,int>())
            .def(py::init<int, int, const EvaporationParameters>())
            .def_readonly("A",&prefragment::A)
            .def_readonly("Z",&prefragment::Z)
            .def_readwrite("atilda",&prefragment::atilda)
            .def_readwrite("config",&prefragment::config)
            .def_readwrite("pairing",&prefragment::pairing);

    py::class_<emission_data>(m,"emission_data")
            .def_readonly("G",&emission_data::G);

    py::class_<emission_results>(m,"emission_results")
            .def_readwrite("g",&emission_results::g)
            .def_readwrite("n",&emission_results::n)
            .def_readwrite("p",&emission_results::p)
            .def_readwrite("d",&emission_results::d)
            .def_readwrite("t",&emission_results::t)
            .def_readwrite("a",&emission_results::a)
            .def_readwrite("he3",&emission_results::he3)
            .def_readwrite("imf",&emission_results::imf)
            .def("dict",[](const emission_results &e){
                py::dict r;
                r["Gg"] = e.g.G;
                r["Gn"] = e.n.G;
                r["Gp"] = (e.p.G);
                r["Gd"] = (e.d.G);
                r["Gt"] = (e.t.G);
                r["Ga"] = (e.a.G);
                r["G3he"] = (e.he3.G);
                r["Gimf"] = (e.imf.G);
                return r;})
            .def("__repr__",[](const emission_results &e){
                py::dict r;
                r["Gg"] = e.g.G;
                r["Gn"] = e.n.G;
                r["Gp"] = (e.p.G);
                r["Gd"] = (e.d.G);
                r["Gt"] = (e.t.G);
                r["Ga"] = (e.a.G);
                r["G3he"] = (e.he3.G);
                r["Gimf"] = (e.imf.G);
                return py::str(r);
            });

    py::class_<EvaporationParameters>(m,"EvaporationParameters")
            .def(py::init<>())
            .def_readwrite("Emax",&EvaporationParameters::Emax)
            .def_readwrite("n_removal_scaling",&EvaporationParameters::n_removal_scaling)
            //.def_readwrite("preset",&EvaporationParameters::preset)
            .def_readwrite("config",&EvaporationParameters::config)
            .def_readwrite("density",&EvaporationParameters::density)
            .def_readwrite("barrier",&EvaporationParameters::barrier)
            .def_readwrite("xn_cs",&EvaporationParameters::xn_cs)
            .def_readwrite("excitation_function",&EvaporationParameters::excitation_function);

    py::enum_<level_density_type>(m,"level_density_type")
        .value("GC_GEM", GC_GEM)
        .value("GC_RIPL", GC_RIPL)
        .value("GC_KTUY05", GC_KTUY05)
        .export_values();

    py::enum_<excitation_function_type>(m,"excitation_function_type")
        .value("GS", GS)
        .value("CUSTOM", CUSTOM)
        .export_values();
        


    m.def("make_model",&make_model,"make Glauber model from python dictionary");
    m.def("make_density",&make_density_from_dict,"make Density from python dictionary");
    m.def("make_evaporation_parameters",&make_evaporation_parameters,"make EvaporationParameters from python dictionary");
    m.def("density_type_to_string",py::overload_cast<const DensityType&>(&density_type_to_string),"DensityType to string");
    m.def("default_nucleus",&defaults_nucleus_python,"get default nucleus from symbol");
    m.def("density_object",&density_object,"get python dictionary of DensityType");
    m.def("nucleus_object",&nucleus_object,"get python dictionary of Nucleus");
    m.def("nucleus_symbol",&nucleus_symbol);
    m.def("sigma_r",py::vectorize(py_sigma_r),"get Sigma_R of the model for energy or array of energies ",py::arg("model"),py::arg("energy"));
    m.def("sigma_cc",py::vectorize(py_sigma_cc),"get Sigma_CC of the model for energy or array of energies",py::arg("model"),py::arg("energy"));
    m.def("r_rms",Rrms<DensityFermi>);
    m.def("r_rms",Rrms<DensityHO>);
    m.def("r_rms",Rrms<DensityType>);
    m.def("fermi_energy",py::vectorize(fermi_energy),"get fermi energy at radius");
    m.def("fermi_momentum",py::vectorize(fermi_momentum),"get fermi momentum at radius");
    //m.def("fermi_momentum_monitz",fermi_momentum_monitz,"get fermi momentum accroding to Monitz et al.");
    m.def("Sp",&ame16::Sp,"Sp");
    m.def("Sn",&ame16::Sn,"Sn");
    m.def("Sa",&ame16::Sa,"Sa");
    m.def("get_nuclear_mass",&ame16::get_nuclear_mass,"get_nuclear_mass");
    m.def("rho_gs",&rho_gs,"rho_gs");
    m.def("w_gs",&w_gs,"w_gs");
    m.def("w_ericson",&w_ericson,"w_ericson");
    m.def("cdf_w_ericson",&cdf_w_ericson,"cdf_w_ericson");
    m.def("cdf_w_gs",&cdf_w_gs,"cdf_w_gs");
    m.def("epax3",&epax::epax3,"epax cross section in mb", py::arg("a_p"), py::arg("z_p"), py::arg("a_t"), py::arg("z_t"),py::arg("a_f"), py::arg("z_f"));
    //m.def("epax_xn_ratios",&epax_xn_ratios,"epax_xn_ratios");
   // m.def("charge_evaporation_probability",&charge_evaporation_probability,"charge evaporation probability",py::arg("A"),py::arg("Z"),py::arg("Ex"),py::arg("h"),py::arg("config")=EvaporationParameters());
    //m.def("charge_evaporation_probability_total",&charge_evaporation_probability_total,"charge evaporation probability including daughters",py::arg("A"),py::arg("Z"),py::arg("Ex"),py::arg("h"),py::arg("config")=EvaporationParameters());
    //m.def("charge_evaporation_probability_simple",&charge_evaporation_probability_simple);
    //m.def("neutron_evaporation_probability",&neutron_evaporation_probability,"neutron evaporation probability",py::arg("A"),py::arg("Z"),py::arg("Ex"),py::arg("h"),py::arg("config")=EvaporationParameters());
    m.def("coulomb_barrier",&coulomb_barrier,"coulomb_barrier");
    m.def("nuclear_potential",&nuclear_potential,"nuclear potential including Coulomnb barrier as a funciton of r, V = Vc-Vn");
    m.def("fusion_barrier",&fusion_barrier,"fusion barrier including Coulomnb barrier as a funciton of r, V = Vc-Vn");
    //m.def("get_evaporation_data",get_evaporation_data);
    //m.def("charge_evaporation_function",&py_charge_evaporation_function);
    m.def("evaporation_ratios",py::overload_cast<prefragment&, double, double>(&evaporation_ratios));
    m.def("evaporation_ratios",py::overload_cast<int, int, double, double, const EvaporationParameters>(&evaporation_ratios),py::arg("A"),py::arg("Z"),py::arg("Ex"),py::arg("j")=0.0,py::arg("evaporation_parameters")=default_evaporation);
    //m.def("level_density",py::overload_cast<prefragment&, double, double>(&level_density),"returns level density and temperature");
    m.def("level_density_gem",&level_density_gem,"returns level density");
    m.def("level_density_kawano",&level_density_kawano,"returns level density");
    m.def("level_density_ripl",&level_density_ripl,"returns level density");
    m.def("Emax",&Emax,"returns default GS excitation function Emax parameter for given nucleus and EvaporationParameters",py::arg("projectile"),py::arg("par"));
    //m.def("width",&width_e);
    //m.def("const_temperature_parameter",py::overload_cast<const prefragment&>(&constant_temperature_parameter));
    //m.def("temperature_parameter",&temperature_parameter);
    //m.def("superfluid_phase_critical_energy",py::overload_cast<const prefragment&, double>(&superfluid_phase_critical_energy));
    //m.def("fermi_gas_density",py::overload_cast<const prefragment&, double>(&fermi_gas_density));
    //m.def("const_temperature_density",py::overload_cast<const prefragment&, double>(&const_temperature_density));
    //m.def("energy_corrected",py::overload_cast<const prefragment&, double>(&energy_corrected));
    //m.def("pairing_energy",&pairing_energy);
    m.def("S",&S);
    m.def("C",&C);
    m.attr("neutron_mass") = neutron_mass;
    m.attr("proton_mass") = proton_mass;
    m.def("custom_excitation_set",&custom_excitation_set,"set excitation function for i hole",py::arg("energies"),py::arg("w"), py::arg("i"));
    m.def("custom_excitation_w",&custom_excitation_w,"get excitation function for given energy and i holes",py::arg("energy"), py::arg("i"));
    m.def("custom_excitation_reset",&custom_excitation_reset,"reset custom excitation function");
    m.def("number_density_cm2",&number_density_cm2,"returns number density per cm2 from thickness per cm2 and molar mass", py::arg("thickness"), py::arg("molar_mass"));
    m.doc() = "Nurex library";
}
