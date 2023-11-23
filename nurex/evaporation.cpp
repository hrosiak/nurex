#include "nurex/evaporation.h"
#include "nurex/epax.h"
#include "nurex/Utils.h"
#include "nurex/data_CGP.h"

using std::pow;

namespace nurex {
EvaporationParameters default_evaporation;
ExcitationFunction custom_excitation;

double rho_gs(double e, int a, const double Emax){
    assert(a>0);
    assert(Emax>0);
    double g0 = 2/Emax;
    //g0 = 16.0;
    double g1 = g0/Emax;
    double sum = 0.0;
    for(int m=0;m<a+1;m++){
        double c = pow(g1,m)/factorial(m);
        c*= pow(g0,a-m)/factorial(a-m);
        c*=pow(e,a+m-1)/factorial(a+m-1);
        assert(c>=0.0);
        if(m%2 != 0) c = -1*c;
        sum += c;
    }
    return sum;
}

double w_gs(double e, int a, const double Emax){
    auto fx = [&](double x)->double{
        return rho_gs(x,a,Emax);
    };

    double hi = Emax;
    if(a>1){
        const double max = Emax*(0.5+a*0.5);
        hi = bisection(fx,Emax,max);
        #ifndef NDEBUG
        if(hi>max || hi<Emax){
            printf("a=%d, hi: %lf, Emax=%lf\n",a,hi, max);
        }
        assert(hi>Emax);
        assert(hi<max);
        #endif
    }
    double norm = integratorGL8.integrate(fx,0,hi);
    if(e>=hi || fx(e)<0.0)return 0.0;
    return fx(e)/norm;
}

double cdf_w_gs(double e, int a, const double Emax){
    auto fx = [&](double x)->double{
        return rho_gs(x,a,Emax);
    };

    double hi = Emax;
    if(a>1){
        hi = Emax*a*0.7;
        hi = bisection(fx,Emax,hi);
        assert(hi>Emax);
    }

    double norm = integratorGL8.integrate(fx,0,hi);
    if(e>=hi)return 0.0;
    double res = integratorGL8.integrate(fx,e,hi)/norm;
    if(res>1.0)res=1.0; // in case of numerical precision is not enough
    return res;
}

double cdf_wfx_gs(Functional& f, double e, int a, const double Emax){
    auto fx = [&](double x)->double{
        return rho_gs(x,a,Emax);
    };
    auto fxx = [&](double x)->double{
        return f.eval(x)*rho_gs(x,a,Emax);
    };

    double hi = Emax;
    if(a>1){
        hi = Emax*a*0.7;
        hi = bisection(fx,Emax,hi);
        assert(hi>Emax);
    }

    double norm = integrator_adaptive.integrate(fx,0,hi);
    if(e>=hi)return 0.0;
    double res = integrator_adaptive.integrate(fxx,e,hi)/norm;
    if(res>1.0)res=1.0; // in case of numerical precision is not enough
    return res;
}

double cdf_wfx_gs(const std::function<double(double)> &f, double e, int a, const double Emax, double hint){
    auto fx = [&](double x)->double{
        return rho_gs(x,a,Emax);
    };
    auto fxx = [&](double x)->double{
        return f(x)*rho_gs(x,a,Emax);
    };

    double hi = Emax;
    if(a>1){
        hi = Emax*a*0.7;
        hi = bisection(fx,Emax,hi);
        assert(hi>Emax);
    }

    double norm = integrator_adaptive.integrate(fx,0,hi);
    if(e>=hi)return 0.0;
    double res;

    if(hint>e){
        res = integrator_adaptive.integrate(fxx,e,hint,1e-3,0,1)/norm;
        res += integrator_adaptive.integrate(fxx,hint,hi,1e-3,0,2)/norm;
    }
    else {
        res = integrator_adaptive.integrate(fxx,e,hi,1e-3,0,4)/norm;
    }
    //res = integrator_adaptive.integrate(fxx,e,hi,1e-3,0,4)/norm;
    if(res>1.0)res=1.0; // in case of numerical precision is not enough
    return res;
}

double cdf_wfx_custom(const std::function<double(double)>&f, double e, int a){
    double res = 0.0;

    auto fxx = [&](double x)->double{
        return f(x)*custom_excitation.w(x,a);
    };
    res = integrator_adaptive.integrate(fxx,e,custom_excitation.emax(a),1e-3,0,4);

    return res;
    }

double cdf_wfx(const std::function<double(double)>&f, double e, int a, const EvaporationParameters &config, const double Emax, double hint){
    if(config.excitation_function == excitation_function_type::GS){
        return cdf_wfx_gs(f,e, a, Emax, hint);
    }
    else if(config.excitation_function == excitation_function_type::CUSTOM){
        return cdf_wfx_custom(f,e, a);
    }    
    else{
        return 0.0;
    }
}

double rho_ericson(double e, int a, const double Emax){
    double g0 = 1/Emax/a;
    double r = pow(e,a-1)*g0/(factorial(int(a))*factorial(int(a-1)));
    return r;
}

double w_ericson(double e, int a, const double Emax){
    auto fx = [&](double x)->double{
        return rho_ericson(x,a,Emax);
    };

    if(e>a*Emax || fx(e)<0.0)return 0.0;
    double norm = integratorGL8.integrate(fx,0,a*Emax);
    return fx(e)/norm;
}

double cdf_w_ericson(double e, int a, const double Emax){
    if(e>a*Emax)return 0.0;
    auto fx = [&](double x)->double{
        return rho_ericson(x,a,Emax);
    };
    double norm = integratorGL8.integrate(fx,0,Emax*a);
    return integratorGL8.integrate(fx,e,Emax*a)/norm;
}

double ExcitationFunction::w(double e, int a) const{
    if(ex.find(a) != ex.end()){
        if(e<ex.at(a).get_min() || e>ex.at(a).get_max())return 0.0;
        else return norm.at(a)*ex.at(a).eval(e);
    }
    else{
        return 0.0;
    }
}

void ExcitationFunction::set(std::vector<double> energy, std::vector<double>rho, int a){
    for (size_t i = 0; i < rho.size(); i++){
        if(rho[i]<0.0)rho[i]=0.0;
    }
    ex.erase(a);
    norm.erase(a);
    ex.emplace(std::make_pair(a,Interpolator(energy, rho)));    
    double res = integrator_adaptive.integrate(ex.at(a),ex.at(a).get_min(),ex.at(a).get_max(),1e-6,0,2);
    if(res<=0.0){
        norm.emplace(a,0.0);
        }
    else {
        norm.emplace(a,1.0/res);
        }
}

double coulomb_potential(double r, int z1 , int z2, double r1, double r2){
    double R = r1+r2;
    double v=0.0;
    double f = 1.44*z1*z2;
    if(r>=R){
        v = f/r;
      }
    else{
        double rr = r/R;
        double rr2 = rr*rr;
        v = f*(3.0 - rr2)/(2.0*R);
      }
    return v;
}

inline double r_param(double A13){
     double R = 1.28*A13 - 0.76 + (0.8/A13);
   // if(A==1)R=R+3.0;
    return R;
    }

/*
double bass(double r, int A1, int Z1, int A2, int Z2){
    double R1 = r_param(A1);
    double R2 = r_param(A2);
    double C1 = R1*(1-((0.9984*0.9984)/(R1*R1)));
    double C2 = R2*(1-((0.9984*0.9984)/(R2*R2)));
    C1  = 1.16*std::pow(A1,1.0/3.0) - 1.39*std::pow(A1,-1.0/3.0);
    C2  = 1.16*std::pow(A2,1.0/3.0) - 1.39*std::pow(A2,-1.0/3.0);
    double s = r-C1-C2;
    double f = C1*C2/(C1+C2);
    return f/((0.033*std::exp(s/3.5)) + (0.007*std::exp(s/0.65)));
}
*/

double bass(double r, double A1_13, double A2_13){
    double R1 = r_param(A1_13);
    double R2 = r_param(A2_13);
    double C1 = R1*(1-(0.98/(R1*R1)));
    double C2 = R2*(1-(0.98/(R2*R2)));

    double s = r-C1-C2;
    double f = C1*C2/(C1+C2);
/*
    if(A1==1){
      s = r - C2;
      f = C2;
    }
    else if (A2==1){
      s = r - C1;
      f = C1;
    }
*/
    return f/((0.033*std::exp(s/3.5)) + (0.007*std::exp(s/0.65)));
}

double nuclear_potential(double r, int A1, int Z1, int A2, int Z2){
    double A1_13 = std::pow(A1,1./3.);
    double A2_13 = std::pow(A2,1./3.);
    double r1 = (A1==1)?0.8:1.24*A1_13;
    double r2 = (A2==1)?0.8:1.24*A2_13;
    return coulomb_potential(r,Z1,Z2,r1,r2)-bass(r,A1_13, A2_13);
};

double fusion_barrier_parametrized(int A1, int Z1, int A2, int Z2){
    const double A1_13 = std::pow(A1,1./3.);
    const double A2_13 = std::pow(A2,1./3.);
    const double eta = (double)Z1*Z2/(A1_13+A2_13);
    double r = -1.01 + (0.93 * eta)  + ((4.53e-4)*eta*eta);
    return r;
};

double fusion_barrier(int A1, int Z1, int A2, int Z2){
    auto fx = [&](double x){
        return nuclear_potential(x, A1, Z1, A2, Z2);
        };
    double rmax = 5*(std::pow(A1,1.0/3.0)+std::pow(A2,1.0/3.0));
    double r = golden_search(fx,0,rmax);
    return fx(r);
}

double asymptotic_level_density_parameter(int A, double bs){
    return 0.073*A + 0.095*bs*std::pow(A,2.0/3.0);
}

prefragment::prefragment(int A, int Z, const EvaporationParameters config):A(A), Z(Z), config(config){
    auto frdmdata = frdm::get_data(A,Z);
    b2 = frdmdata.b2;
    esp = frdmdata.EFLmic;
    pairing = pairing_energy(A, Z, config.density);

    // asymptotic level density
    if(config.density == level_density_type::GC_KTUY05){
        atilda = (0.1143*A) + (7.6e-5*A*A);
    }
    else if(config.density == level_density_type::GC_GEM || config.density == level_density_type::GC_RIPL){
        if(Z<9 || (A-Z)<9){
            atilda = A/8.0;
        }
        else {
        double aprime = (
        (Z >=54 && Z<78)
        || (Z >=86 && Z<98)
        || ((A-Z) >=86 && (A-Z)<122)
        || ((A-Z) >=130 && (A-Z)<150)
        )? 0.12 : 0.142;
        atilda = A*(aprime + 0.00917*esp);
        }
    }
    else{
        atilda = asymptotic_level_density_parameter(A,bs_ratio(b2,frdmdata.b4));
    }

}

double S(prefragment &f, int Ap, int Zp){
    if(Ap == 1 && Zp==0){
        if(f.Sn<=0.)f.Sn = ame16::Sn(f.A,f.Z);
        return f.Sn;
    }
    else if(Ap == 1 && Zp ==1){
        if(f.Sp<=0.)f.Sp = ame16::Sp(f.A,f.Z);
        return f.Sp;
    }
    else if(Ap == 2 && Zp ==1){
        if(f.Sd<=0.)f.Sd = ame16::S(f.A,f.Z,1,1);
        return f.Sd;
    }
    else if(Ap == 3 && Zp ==1){
        if(f.St<=0.)f.St = ame16::S(f.A,f.Z,1,2);
        return f.St;
    }
    else if(Ap == 3 && Zp ==2){
        if(f.She3<=0.)f.She3 = ame16::S(f.A,f.Z,2,1);
        return f.She3;
    }
    else if(Ap == 4 && Zp ==2){
        if(f.Sa<=0.)f.Sa = ame16::Sa(f.A,f.Z);
        return f.Sa;
    }
    else{
        return ame16::S(f.A, f.Z, Zp, Ap-Zp);
    }
}

double C(prefragment &f, int Ap, int Zp){
    if(f.config.barrier == barrier_type::none)return 0.0;
    if(Ap == 1 && Zp==0)return 0.;
    else if(Ap == 1 && Zp ==1){
        if(f.Cp<=0.0){
            f.Cp = (f.config.barrier == barrier_type::parametrized)?
                       fusion_barrier_parametrized(f.A-1,f.Z-1,1,1)
                       :fusion_barrier(f.A-1,f.Z-1,1,1);
        }
        return f.Cp;
    }
    else if(Ap == 2 && Zp ==1){
        if(f.Cd<=0.0){
            f.Cd = (f.config.barrier == barrier_type::parametrized)?
                        fusion_barrier_parametrized(f.A-2,f.Z-1,2,1)
                      :fusion_barrier(f.A-2,f.Z-1,2,1);
        }
        return f.Cd;
    }
    else if(Ap == 3 && Zp ==1){
        if(f.Ct<=0.0){
            f.Ct = (f.config.barrier == barrier_type::parametrized)?
                       fusion_barrier_parametrized(f.A-3,f.Z-1,3,1)
                      :fusion_barrier(f.A-3,f.Z-1,3,1);
        }
        return f.Ct;
    }
    else if(Ap == 3 && Zp ==2){
        if(f.Che3<=0.0){
            f.Che3 = (f.config.barrier == barrier_type::parametrized)?
                     fusion_barrier_parametrized(f.A-3,f.Z-2,3,2)
                      :fusion_barrier(f.A-3,f.Z-2,3,2);
        }
        return f.Che3;
    }
    else if(Ap == 4 && Zp ==2){
        if(f.Ca<=0.0){
            f.Ca = (f.config.barrier == barrier_type::parametrized)?
                   fusion_barrier_parametrized(f.A-4,f.Z-2,4,2)
                   :fusion_barrier(f.A-4,f.Z-2,4,2);
        }
        return f.Ca;
    }
    else{
        return  (f.config.barrier == barrier_type::parametrized)?
                    fusion_barrier_parametrized(f.A-Ap, f.Z-Zp,Ap, Zp)
                   :fusion_barrier(f.A-Ap, f.Z-Zp,Ap, Zp);
    }
}

double charge_evaporation_probability_simple(int A, int Z, double Ex, int a){
    assert(a>0);
    assert(Ex>0);
    assert(A-a>0);
    assert(A>4);
    double Sp = ame16::Sp(A, Z);
    double Sa = ame16::Sa(A, Z);
    double Cp = fusion_barrier(A-1, Z-1, 1, 1);
    double Ca = (A-4>0)?fusion_barrier(A-4, Z-2, 4, 2):99999999;
    double V = std::min(Sp + Cp, Sa + Ca);
    if(V<=0.1)return 1.0;
    return cdf_w_gs(V, a, Ex);
}

double charge_evaporation_probability_total(int A, int Z, double Emax, int a, const EvaporationParameters &config){
    assert(a>0);
    assert(A-a>0);
    assert(A>4);
    prefragment CN(A,Z, config);
    double Sn = S(CN,1,0);
    double Sp = S(CN,1,1);
    double Sa = S(CN,4,2);
    double Cp = C(CN,1,1);
    double Ca = (A-4>0)?C(CN,4,2):99999999;
    double V = std::min(Sp + Cp, Sa + Ca);
    if(V<=0.1)return 1.0;
    if(Emax<=0.)return 0.0;

    auto jd = angular_momentum_distribution(A+a,A,0.0);
    auto frch = [&](double x)->double{
        //return charge_evaporation_function(CN, x,static_cast<int>(jd.first));
        return charge_evaporation_function(CN, x, 0);
    };
    //return cdf_wfx_gs(frch,V, a, Emax,(Sn>V)?Sn:0.0);
    return cdf_wfx(frch,V, a, config, Emax,(Sn>V)?Sn:0.0);
}

emission_data get_emission_data(prefragment &f, int Ap, int Zp, double Ex, double jd){
    const double B = C(f,Ap,Zp);
    const double SE = S(f,Ap,Zp);
    const double Eeff = Ex-SE-B;
    emission_data r;
    double w;
    if(Eeff<0.0)return r;
    int Af = f.A - Ap;
    int Zf = f.Z - Zp;
    prefragment frag(Af,Zf, f.config);
    double Erot = (jd>0.0)?0.5*jd*jd/J(Af, frag.b2).first:0.0;
    w = width_e(f,Ap,Zp,B, B+SE,Ex-Erot,jd);
    r.G = w;
    //r.rho = g_d.first;
    //r.T = g_d.second;
    return r;
}

emission_data mean_emission_data(prefragment &f, int Ap, int Zp, double Ex, double lm){
        const double B = C(f,Ap,Zp);
        const double SE = S(f,Ap,Zp);
        const double Eeff = Ex-SE-B;
        double G = 0.;
        double rho;
        double T;
        if(Eeff<0.1)return {};
        if(lm>0.5){
            auto vj = get_discrete_vector_by_area( l_orb_distribution(f.A, f.A-Ap, lm, Eeff, f.atilda) );
            for (auto& entry : vj){
               int jm = static_cast<int>(lm);
               auto data = get_emission_data(f,Ap, Zp, Ex,std::abs(jm-entry.first));
               G+=entry.second*data.G;
               T+=entry.second*data.T;
               rho+=entry.second*data.rho;
            }
        }
        else{
            auto data = get_emission_data(f,Ap, Zp, Ex,0.0);
            G = data.G;
            T = data.T;
            rho = data.rho;
        }
        return {G,rho,T};
}

double charge_evaporation_function(prefragment &f, double Ex, double j){
    const double T_f = T_freezeout(f.A);
    const double T_in = std::sqrt(Ex/f.atilda);
    if(T_in-T_f > 0.1)return 1.0; // break-up
    const double Sp = S(f,1,1);
    const double Sa = S(f,4,2);
    if(Sp<0. || Sa<0.)return 1.0;
    const double V = std::min(Sp + C(f,1,1), Sa + C(f,4,2));

    if(V<=0.1)return 1.0;   // charge barrier is low
    if(Ex<V){               // Excitation function below charge barrier so no charge evaporation
        return 0.0;
    }
    else{
//        if(Ex<Sn)return 1.0; // no neutron evap.
        //auto g_CN = level_density(f, Ex, j).first;
        auto data = evaporation_ratios(f,Ex,j);
        auto Gn = data.n;
        auto Gp = data.p;
        auto Gd = data.d;
        auto Gt = data.t;
        auto G3he = data.he3;
        auto Ga = data.a;
        emission_data Gimf = data.imf;
        double sum = (Gn.G+Gp.G+Ga.G+Gd.G+Gt.G+G3he.G+Gimf.G);
        if(sum==0.0)return 0.0;

        prefragment fn(f.A-1, f.Z, f.config);
        double Sn = S(f,1,0);
        double Vf = std::min(S(fn,1,1) + C(fn,1,1), S(fn,4,2) + C(fn,4,2));
        double P = 1.0;
        if(Sn>0.1){
            double T = level_density(f, Ex, j).second;
            //double T = Gn.T; // temperature of neutron daughter
            Ex -= Sn + 2*T*sqrt(2./nurex::PI);
            if( Ex < Vf )P = 0.0;
            else P = charge_evaporation_function(fn,Ex,0); //j set to 0 ?
        }
        return 1.- ((1.-P)*Gn.G/sum);
    }
}

double neutron_evaporation_probability(int A, int Z, double Ex, int a, const EvaporationParameters &config){
    assert(a>0);
    assert(Ex>0);
    assert(A-a>0);
    assert(A>4);
    prefragment CN(A,Z, config);
    double Sn = S(CN,1,0);
    if(Sn<=0.1)return 1.0;
    if(Ex<Sn)return 0.0;    
    auto frn = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,0);

        prefragment fn(A-1, Z, config);
        double P = 1.0;
        double V = std::min({S(fn,1,1) + C(fn,1,1), S(fn,4,2) + C(fn,4,2), S(fn, 1, 0)});
        double T = level_density(CN, Ex, 0).second;
        double Sn = S(CN,1,0);
        double ex = x - Sn - 2*T*sqrt(2./nurex::PI);
        
        if(V>0.1){
            if( ex < V )
                P = 0.0;    
            else{                
                auto data2 = evaporation_ratios(fn,ex,0);
                P = 1.0 - data2.g.G;
            } 
        }     
        if(!test) P=0.0;
        assert(P>=0.0 && P<=1.0);
        //if(P>0.0)printf("%lf %lf %lf\n",P, V, ex);
        return (1.0 - P)*data.n.G;
    };

    //return cdf_wfx_gs(frn,Sn, a, Ex);
    return cdf_wfx(frn,Sn, a, config, Ex);
}

double charge_evaporation_probability(int A, int Z, double Emax, int a, const EvaporationParameters &config){
    assert(a>0);
    assert(Emax>0);
    assert(A-a>0);
    assert(A>4);
    prefragment CN(A,Z, config);
    const double V = std::min(S(CN,1,1) + C(CN,1,1), S(CN,4,2) + C(CN,4,2));
    if(V<=0.1)return 1.0;
    if(Emax<V)return 0.0;
    auto jd = angular_momentum_distribution(A+a,A,0.0);
    auto frn = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return 1.0-data.n.G-data.g.G;
    };

    //return cdf_wfx_gs(frn,V, a, Emax);
    return cdf_wfx(frn,V, a, config, Emax);
}

double total_evaporation_probability(int A, int Z, double Ex, int a, const EvaporationParameters &config){
    assert(a>0);
    assert(Ex>0);
    assert(A-a>0);
    assert(A>4);
    prefragment CN(A,Z, config);
    double V = std::min(S(CN,1,1) + C(CN,1,1), S(CN,4,2) + C(CN,4,2));
    V = std::min(S(CN,1,0), V);
    if(V<=0.1)return 1.0;
    if(Ex<V)return 0.0;
    auto jd = angular_momentum_distribution(A+a,A,0.0);
    auto frn = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return 1.0-data.g.G;
    };

    //return cdf_wfx_gs(frn,V, a, Ex);
    return cdf_wfx(frn,V, a, config, Ex);
}

emission_results evaporation_probabilities(int A, int Z, double Ex, int a, const EvaporationParameters &config){
    assert(a>0);
    assert(Ex>0);
    assert(A-a>0);
    assert(A>4);
    emission_results res; 
    prefragment CN(A,Z, config);
    double V = std::min(S(CN,1,1) + C(CN,1,1), S(CN,4,2) + C(CN,4,2));
    V = std::min(S(CN,1,0), V);
    if(Ex<V){
        res.g.G = 1.0;
        return res;
    };
    auto jd = angular_momentum_distribution(A+a,A,0.0);
    
    auto frn = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return data.n.G;
    };

    auto frp = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return data.p.G;
    };

    auto frd = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return data.d.G;
    };

    auto frt = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return data.t.G;
    };

    auto fra = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return data.a.G;
    };

    auto frhe3 = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return data.he3.G;
    };

    auto frimf = [&](double x)->double{
        auto data = evaporation_ratios(CN,x,static_cast<int>(jd.first));
        return data.imf.G;
    };
    

    res.n.G =  cdf_wfx(frn,V, a, config, Ex);
    res.p.G =  cdf_wfx(frp,V, a, config, Ex);
    res.d.G =  cdf_wfx(frd,V, a, config, Ex);
    res.t.G =  cdf_wfx(frt,V, a, config, Ex);
    res.he3.G =  cdf_wfx(frhe3,V, a, config, Ex);
    res.a.G =  cdf_wfx(fra,V, a, config, Ex);
    res.imf.G =  cdf_wfx(frimf,V, a, config, Ex);
/*
    if(config.excitation_function == excitation_function_type::GS){
        res.n.G =  cdf_wfx_gs(frn,V, a, Ex);
        res.p.G =  cdf_wfx_gs(frp,V, a, Ex);
        res.d.G =  cdf_wfx_gs(frd,V, a, Ex);
        res.t.G =  cdf_wfx_gs(frt,V, a, Ex);
        res.he3.G =  cdf_wfx_gs(frhe3,V, a, Ex);
        res.a.G =  cdf_wfx_gs(fra,V, a, Ex);
        res.imf.G =  cdf_wfx_gs(frimf,V, a, Ex);
    }
    else if(config.excitation_function == excitation_function_type::CUSTOM){
        res.n.G =  cdf_wfx_custom(frn,V, a);
        res.p.G =  cdf_wfx_custom(frp,V, a);
        res.d.G =  cdf_wfx_custom(frd,V, a);
        res.t.G =  cdf_wfx_custom(frt,V, a);
        res.he3.G =  cdf_wfx_custom(frhe3,V, a);
        res.a.G =  cdf_wfx_custom(fra,V, a);
        res.imf.G =  cdf_wfx_custom(frimf,V, a);
    }
  */  
    return res;
}

emission_results evaporation_ratios(prefragment &f, double Ex, double j){
    emission_results r;
    const double T_f = T_freezeout(f.A);
    const double T_in = std::sqrt(Ex/f.atilda);
    //double Efreeze = T_f*T_f*atilda;

    if(T_in-T_f > 0.1){ // break-up
        r.imf.G = 1.0;
        return r;
    }

    const double Sn = S(f,1,0);
    const double V = std::min(S(f,1,1) + C(f,1,1), S(f,4,2) + C(f,4,2));

    if(Ex<V && Ex<Sn){
        r.g.G = 1.0;
        return r;} // energy under particle emission threshold
    else if(Ex<V && Ex>Sn){
        //auto Gn = mean_emission_data(f,1,0,Ex,j);
        r.n.G=1.0;
        return r;
    }   // only neutron emission possible?
    else{
        //auto g_CN = level_density(f, Ex, j).first;

        auto Gn = ( (f.config.config&evaporation_config_type::disable_neutron)==0 )?mean_emission_data(f,1,0,Ex,j):emission_data{};
        auto Gp = mean_emission_data(f,1,1,Ex,j);
        auto Gd = mean_emission_data(f,2,1,Ex,j);
        auto Gt = mean_emission_data(f,3,1,Ex,j);
        auto G3he = mean_emission_data(f,3,2,Ex,j);
        auto Ga = mean_emission_data(f,4,2,Ex,j);
        emission_data Gimf;

        if( (f.config.config&evaporation_config_type::disable_imf)==0 ){
        for(int zz=3;zz<=f.Z/2;zz++){
            double Gz = 0.0;
            for (int nn=3;nn<f.A-zz-1;nn++) {
                int aa = zz+nn;
                if(f.A <=2*aa)continue;
                if(ame16::Sp(aa,zz)<=0. || ame16::Sp(f.A-aa,f.Z-zz)<=0.)continue; // check if residual and evaporated are bound
                if(ame16::Sn(aa,zz)<=0. || ame16::Sn(f.A-aa,f.Z-zz)<=0.)continue; // if unbound its not considered
                Gz += mean_emission_data(f,aa,zz,Ex,j).G;
            }
            if(Gz==0.0)break;
            Gimf.G += Gz;
        }
        }
        else{
            Gimf.G = 0.0;
        }

        double sum = (Gn.G+Gp.G+Ga.G+Gd.G+Gt.G+G3he.G+Gimf.G);
        if(sum==0.0){
            r.g.G = 1.0;
            return r;
            } //no calculated emmission possible, put is as gamma

        Gn.G = Gn.G/sum;
        Gp.G = Gp.G/sum;
        Gd.G = Gd.G/sum;
        Gt.G = Gt.G/sum;
        G3he.G = G3he.G/sum;
        Ga.G = Ga.G/sum;
        Gimf.G = Gimf.G/sum;
        assert(Gn.G < 1.0001);assert(Gp.G < 1.0001);assert(Gd.G < 1.0001);
        assert(Gt.G < 1.0001);assert(G3he.G < 1.0001);assert(Ga.G < 1.0001);
        assert(Gimf.G < 1.0001);


        r.g.G = 0.0;
        r.n = Gn;
        r.p = Gp;
        r.d = Gd;
        r.t = Gt;
        r.he3 = G3he;
        r.a = Ga;
        r.imf = Gimf;

        return r;
    }
}



std::pair<double, double> J(int A, double b2){
    const double C = 0.4*std::pow(A,5.0/3.0)*1.16*1.16*nurex::atomic_mass_unit/(nurex::hc_transp*nurex::hc_transp);
    double Jper = C*(1.0+0.50*b2*sqrt(5.0/(4.0*nurex::PI)));
    double Jpar = C*(1.0-b2*sqrt(5.0/(4.0*nurex::PI)));
    return {Jper, Jpar};
}

double bs_ratio(double b2, double b4){
    double bs = 1.0;
    double a2 = a2_from_b2(b2);
    double a4 = a4_from_b4(b4);
    bs = 1.0 + 0.4*a2*a2 - 4.0/105.0*a2*a2*a2 - 66.0/175.0*a2*a2*a2*a2 - 4.0/35.0*a2*a2*a4 + a4*a4; //6.22
    return bs;
}

double shell_effect_damping(int A, double a, double Ex){
    const double gamma =   2.5 * a * std::pow(A,-4.0/3.0);
    return  1-std::exp(-gamma*Ex);
}

double superfluid_phase_critical_energy(const prefragment &f, double Ex){
    const double gamma =   2.5 * f.atilda * std::pow(f.A,-4.0/3.0);
    const double dU = f.esp;
    const double C= sqrt(1.0+gamma*dU);
    const double P = 17.60*(pow(f.A,-0.699));
    return f.atilda*P*C*P*C;
    //return 10.0;
}

double constant_temperature_parameter(const prefragment &f){
    if(f.config.density==level_density_type::GC_GEM){
        double Ux = 2.5 +150.0/f.A;
        double inv = std::sqrt(f.atilda/Ux) - 1.5/Ux;
        return 1.0/inv;
    }
    else{
        double gamma =   2.5 * f.atilda * std::pow(f.A,-4.0/3.0);
        return 17.60*std::pow(f.A,-0.699) * std::sqrt(1.+gamma*f.esp);
    }
}

double const_temperature_density(const prefragment &f, double Ex){
    if(f.config.density==level_density_type::GC_GEM){
        double Ux = 2.5 +150.0/f.A;
        double T = 1.0/(std::sqrt(f.atilda/Ux) - 1.5/Ux);
        double Ecrit = Ux + pairing_energy(f.A,f.Z,1);
        double E0 = Ecrit - T * (log(T) - 0.25*log(f.atilda) - 1.25*log(Ux) + 2*sqrt(f.atilda*Ux));
        return (PI/12.0)*exp( (Ex - E0)/T)/T;
    }
    else if(f.config.density==level_density_type::GC_RIPL){
        double gamma = 0.325*pow(f.A,-1.0/3.0);
        const double T = 17.6*pow(f.A,-0.699)*sqrt(1.0 + gamma*f.esp);
        double U0 = -0.079*f.esp;
        int N = f.A - f.Z;
        int P = 0;
        if(f.Z%2 == 0)P++;
        if(N%2 == 0)P++;
        if(P == 0)U0 += -11.17*std::pow(f.A, -0.464) + 0.285;
        else if(P == 2)U0 += 11.17*std::pow(f.A, -0.464) - 0.520;
        else U0 += -0.39 - 0.00058*f.A;
        return  exp( (Ex - U0)/T)/T;
    }
    else{
        double Tct = constant_temperature_parameter(f);
        double E = std::max(Ex, 0.0);
        double U0 = -0.079*f.esp;
        int N = f.A - f.Z;
        int P = 0;
        if(f.Z%2 == 0)P++;
        if(N%2 == 0)P++;
        if(P == 0)U0 += -11.17*std::pow(f.A, -0.464) + 0.285;
        else if(P == 2)U0 += 11.17*std::pow(f.A, -0.464) - 0.520;
        else U0 += -0.39 - 0.00058*f.A;
        return std::exp((E - U0)/Tct)/Tct;
    }
    }

double pairing_energy(int A, int Z, int conf){
    const int N = A - Z;
    int P = 0;
    int res = 0;
    if(Z%2 == 0)P++;
    if(N%2 == 0)P++;


    if(conf == level_density_type::GC_GEM || conf == level_density_type::GC_RIPL){
        return cgc_pairing_energy(A,Z);
    }
    else if(conf == 10){
        double M1 = ame16::get_nuclear_mass(A+2,Z+1);
        double M2 = ame16::get_nuclear_mass(A,Z);
        double M3 = ame16::get_nuclear_mass(A-2,Z-1);
        if(M1>0 && M2>0 && M3>0){
            return 0.25*(M1 - (2*M2) + M3)*atomic_mass_unit;
        }
        else {
            return P*12*std::pow(A, -0.5);
        }
    }
    else if(conf == level_density_type::GC_KTUY05){
        const double c1n = 0.45;
        const double c2n = 1.6;
        const double c1p = c1n-0.05;
        const double c2p = c2n-0.02;
        const double c00 = 0.1070;
        const double c001 = 20.14;
        const double alpha = 0.25;
        const double Aad = 2.0;
        const double rad = 0.1268;
        const double rI = 3.0;
        const double rq = 1.4996;
        const double r0 = 1.04;

        double dn, dp;
        const double C1 = std::pow(3*PI*PI,1./3.);
        const double C2 = std::pow(4.0*PI/3.0, 2.0/3.0);
        double nex = (double)(N-Z)/(double)A;
        double Rn = (r0*std::pow(Aad + A,1.0/3.0)*(1+nex*nex))
                   +rad+ (rI*nex) + (rq*nex*nex);

        dn = 0.5*PI*PI*hc_transp*hc_transp/(neutron_mass*C1*std::pow(N,1.0/3.0)*C2*Rn*Rn);
        dp = 0.5*PI*PI*hc_transp*hc_transp/(proton_mass*C1*std::pow(Z,1.0/3.0)*C2*Rn*Rn);

        double Mn = (N%2)?0.0:c1n*dn + c2n*std::pow(dn,alpha);
        double Mp = (Z%2)?0.0:c1p*dp + c2p*std::pow(dp,alpha);
        double Moo = 0.0;
        if( (N%2)==0 && (Z%2)==0){
            Moo = (N==Z)?0.0:c00*(1+c001*(pow(A,-1.0/3.0) - pow(A,-2.0/3.0)))*Mp*Mn/(Mp+Mn);
        }
        return Mn + Mp - Moo;
    }
    else if(conf == level_density_type::ABLA) {
        if(P == 2) return 22.34*std::pow(A,-0.464)-0.235;
        else if(P == 1) return (0.285+11.17*std::pow(A,-0.464)-0.390-0.00058*A);
        else return 0.0;
    }
    else{
        return P*12*std::pow(A, -0.5);
    }
}

double energy_corrected(int A, int Z, double Ex){
    return energy_corrected(prefragment(A,Z), Ex);
};

double energy_corrected(const prefragment &f, double Ex){
    const double dU = f.esp;
    const double k_eff = shell_effect_damping(f.A, f.atilda, Ex);
    //const double dP = -0.25*144 * 6.0 * f.atilda /(PI*PI*f.A) + 22.34*pow(f.A,-0.464)-0.235;
    const double dP = -0.25*144 * 6.0 * f.atilda /(PI*PI*f.A) + 24.0/std::sqrt(f.A);
    //const double Ecrit = 10.0;
    const double Ecrit = superfluid_phase_critical_energy(f, Ex);
    const double heff = (Ex>Ecrit)?1.0:1.0-std::pow(1-(Ex/Ecrit),2);
    return Ex + (dU*k_eff) + (dP*heff);
};

double fermi_gas_density(int A, int Z, double Ex){
    return fermi_gas_density(prefragment(A,Z), Ex);
}

double fermi_gas_density(const prefragment &f, double Ex){
    double Ecor = energy_corrected(f, Ex);
    double a = f.atilda * Ecor/Ex;

    //double d = std::pow(f.atilda,0.25) * std::pow(Ex,5./4.);
    //double S = 2*std::sqrt(f.atilda*Ecor);
    double d = std::pow(a,0.25) * std::pow(Ex,5./4.);
    double S = 2*std::sqrt(a*Ecor);
    return 0.1477045*std::exp(S)/d;
}

double level_density_kawano(prefragment &f, double E, double Em){
    const double Delta = f.pairing;
    const double Eeff = E - Delta;
    double k = shell_effect_damping(f.A, f.atilda,Eeff);
    double a = f.atilda*(1+(f.esp/Eeff)*k);

    const double Ux = 2.5 + 150.0/f.A;
    if(Em<=0.0)Em = Ux + Delta;

    double res = 0;
    double resfg=0.0;
    if(E>=Em && Eeff>0.5 && a>0){
        resfg = (1.0/12.0)*exp(2.0*sqrt(a*Eeff))/(std::pow(a,0.25)*std::pow(Eeff,5.0/4.0));
    }
    if(E>30)res = resfg;
    else {
        double T = 48.07*std::pow(f.A,-0.8834)*sqrt(1-(0.1*f.esp));
        double fA;
        if(closest_shell_difference(f.A) || closest_shell_difference(f.Z))fA = -0.007473*f.A - 1.725;
        else if( (f.A%2)==0 && (f.Z%2)==0)fA = -0.01525*f.A - 1.190;
        else if( ((f.A-f.Z)%2) >0 && (f.Z%2)>0)fA = -0.01416*f.A - 1.595;
        else fA = -0.01525*f.A - 1.190;
        double E0 = Delta - (0.16*f.esp) + (T*fA);
        res = exp( (E - E0)/T)/T;
        if(resfg && resfg<res)res = resfg;
    }
    return res;
}

double level_density_gem(prefragment &f, double E){
    const double Delta = f.pairing;
    const double Ux = 2.5 + 150.0/f.A;
    const double Ex = Ux + Delta;
    double a = f.atilda;
    double u = 0.05*(E-Delta);
    double Ai = f.A * (0.1375 - 8.36e-5*f.A);
    double c = (1.0-exp(-u))/u;
    a = (f.atilda*c) + (Ai*(1-c));

    double res = 0;
    double Eeff = E-Delta;
    if(E>=Ex && Eeff>0.0 && a>0.0){
        res = exp(2.0*sqrt(a*Eeff))/(std::pow(a,0.25)*std::pow(Eeff,5.0/4.0));
    }
    else {
        const double T = 1.0/(std::sqrt(f.atilda/Ux) - 1.5/Ux);
        const double E0 = Ex - T * (log(T) - 0.25*log(a) - 1.25*log(Ux) + 2*sqrt(a*Ux));
        res = exp( (E - E0)/T)/T;
    }
    return (PI/12.0)*res;
}

double level_density_ripl(prefragment &f, double E, double Em){
    const double Delta = f.pairing;
    const double Ux = 2.5 + 150.0/f.A;
    const double Ex = (Em>0.0)?Em:Ux + Delta;
    double a = f.atilda;
    f.atilda = 0.0959*f.A + 0.1468*pow(f.A, 2./3.);
    double gamma = 0.325*pow(f.A,-1.0/3.0);
    a = f.atilda*(1+f.esp*(1.0-exp(-gamma*E))/E);
    if(a<=0.0)return 0.0;
    double res = 0, resfg = 0.0;
    double Eeff = E-Delta;


    if(E>=Ex && Eeff>0.0 && a>0.0){
        resfg = (sqrt(PI)/12.0)*exp(2.0*sqrt(a*Eeff))/(std::pow(a,0.25)*std::pow(Eeff,5.0/4.0));
    }
    if(E>30){
        res = resfg;
    }
    else {
        const double T = 17.6*pow(f.A,-0.699)*sqrt(1.0 + gamma*f.esp);
        double U0 = -0.079*f.esp;
        int N = f.A - f.Z;
        int P = 0;
        if(f.Z%2 == 0)P++;
        if(N%2 == 0)P++;
        if(P == 0)U0 += -11.17*std::pow(f.A, -0.464) + 0.285;
        else if(P == 2)U0 += 11.17*std::pow(f.A, -0.464) - 0.520;
        else U0 += -0.39 - 0.00058*f.A;
        res = exp( (E - U0)/T)/T;
        if(resfg>0.0 && resfg<res)res=resfg;
    }
    return res;
}

std::pair<double, double> level_density(prefragment &f, double Ex, double j){
    double rho=0.0;
    double T=0.0;
    const double Ecrit = superfluid_phase_critical_energy(f, Ex);
    const double Eeff = Ex - f.pairing;
    double Ect = Eeff;
    if(Ex<Ecrit){
        rho = const_temperature_density(f, Ect);
        T = constant_temperature_parameter(f);
    }
    else if(Ex>30){
        rho = fermi_gas_density(f, Eeff);
        T = temperature_parameter(Eeff,f.atilda);
    }
    else{
        rho = const_temperature_density(f, Ect);
        T = constant_temperature_parameter(f);
        double tfg = temperature_parameter(Eeff, f.atilda);
        double rhofg = fermi_gas_density(f, Eeff);
        if(rhofg<rho){rho = rhofg;}
        if(tfg>T){T = tfg;}
    }
    double Jper = 1.0;
    double Jpar = 1.0;
    double sigma2_per = 1.0;
    double sigma2 = 1.0;
    double k = 1.0;
    double j_factor = 1.0;

    if(   (j>0 && f.Z>56 && (evaporation_config_type::disable_k_factor&f.config.config) == 0)
        ||(f.config.config&evaporation_config_type::j_sum)
        ||(f.config.config&evaporation_config_type::j_single)){
        double Jper = J(f.A, f.b2).first;
        double Jpar = J(f.A, f.b2).second;
        double sigma2_per = Jper*T;
        double sigma2 = sigma2_per + (Jpar*T);
    }

    if(j>0 && f.Z>56 && (evaporation_config_type::disable_k_factor&f.config.config) == 0){
        double Erot = 0.5*j*j/std::sqrt((Jper*Jper) + (Jpar*Jpar) );
        //k = Krot(f.A,f.Z, Ecor-Erot, sigma2_per, f.b2);
        k = Krot(f.A,f.Z, Eeff-Erot, sigma2_per, f.b2);
    }

    if( (f.config.config&evaporation_config_type::j_sum)){
        j_factor = 1./(sqrt(2.*nurex::PI*sigma2));
        }
    else if( (f.config.config&evaporation_config_type::j_single)){
        j_factor = J_density_factor(j, sigma2);
        }
    assert(rho>0);
    assert(T>0);
    return {k*j_factor*rho, T};
};

std::pair<double,double> angular_momentum_distribution(int A, int Af, double beta){
    double mean = 0.16*pow(A,2./3.)*(1.-(2.*beta/3.));
    double sigma = mean*(Af*(A-Af))/(double)(A-1);
    return {sqrt(mean),sqrt(sigma)};
}

std::pair<double, double> l_orb_distribution(int Am, int Ad, double l_m, double Ef, double a){
    double l_orb = -1.0;
    double sl_orb = 0.0;
    assert(Ef>0.);
    if(l_m <0.1)return{0.0,0.0};
    if(Ef <0.1)return{-1.0,0.0};
    const int Af = Am - Ad;
    const double sqrts4 = sqrt(a*Ef);
    const double C = atomic_mass_unit*1.3*1.3/(hc_transp*c);
    const double theta_m = 0.4* C * std::pow(Am,5./3.);
    //const double theta_orb = C*std::pow(std::pow(Af,1./3.)+std::pow(Ad,1./3.),2)*( (Af*Af*Ad) + (Ad*Ad*Af))/(Am*Am);
    const double theta_orb = 0.0323*std::pow(std::pow(Af,1./3.)+std::pow(Ad,1./3.),2)*( (Af*Af*Ad) + (Ad*Ad*Af))/(Am*Am);
    l_orb = theta_orb * ( (l_m / theta_m) + (sqrts4 /(a*l_m)));
    sl_orb = sqrt(sqrts4*theta_orb/a);
    return {l_orb, sl_orb};
}
/*
std::pair<double, double> l_orb(int Am, int Ad, double l_m, double Ef, double a){
    auto par = l_orb_distribution(Am, Ad,l_m, Ef, a);
    auto fx = [](double x){
        return
    }
};
*/

double Krot(int A, int Z, double Ex, double sigma2, double b2){
    double Kr = 1.0;
    if(std::abs(b2)<0.15){
      int dn = std::abs(closest_shell_difference(A-Z));
      int dz = std::abs(closest_shell_difference(Z));
      double beff = 0.022 + (0.003*dn) + (0.005*dz);
      sigma2 = 75.0*beff*beff*sigma2;
    }

    if(sigma2>1.0){
        Kr = ((sigma2-1.0)/(1.0+std::exp((Ex - 40.0)/10.0))) + 1.0;
    }

    return Kr;
};

double J_density_factor(double j, double sigma2){
      //return (2.*j+1.)*std::exp(-1.*j*(j+1.0)/(2.0*sigma2))/(std::sqrt(8.0*nurex::PI)*std::pow(sigma2,1.5));
      return (2.*j+1.)*std::exp(-1.*(j+0.5)*(j+0.5)/(2.0*sigma2))/(std::sqrt(8.0*nurex::PI)*std::pow(sigma2,1.5));
}

double sigma_c(int Ad, int Ap, double epsilon, double Vb){
    double Rg = 1.16 * (std::pow((double)Ad,1./3.) + std::pow((double)Ap,1./3.));
    double mu = (double)(Ad*Ap)/(Ad+Ap);
    double Ecm = epsilon * (Ad-Ap)/(double)Ad;
    double Ra = hc_transp * std::sqrt(0.5/(mu*Ecm));
    return PI * (Ra+Rg)*(Ra+Rg)*(1.-Vb/epsilon);
}


double width_gem(int Am, int Zm, int Ap, int Zp, double B, double SB, double Ex, double j){
    double res = 0;
    if(SB>Ex)return 0.0;
    const int Af = Am - Ap;
    const int Zf = Zm - Zp;

    double Rj = (Zp>1)?1.2:0;
    double k;
    if(Zp == 4){
        if(Zf<=20)k=0.81;
        else if(Zf<=30){int dZ = Zf-20; k=0.81 + dZ*0.1*(0.85-0.81);}
        else if(Zf<=40){int dZ = Zf-30; k=0.85 + dZ*0.1*(0.89-0.85);}
        else if(Zf<=50){int dZ = Zf-40; k=0.89 + dZ*0.1*(0.93-0.89);}
        else k=0.93;
    }
    else{
        if(Zf<=20)k=0.51;
        else if(Zf<=30){int dZ = Zf-20; k=0.51 + dZ*0.1*(0.60-0.51);}
        else if(Zf<=40){int dZ = Zf-30; k=0.6 + dZ*0.1*(0.66-0.6);}
        else if(Zf<=50){int dZ = Zf-40; k=0.66 +  dZ*0.1*(0.68-0.66);}
        else k=0.68;
    }

    double Ccoul = k*1.44*Zp*Zf/(1.7*pow(Af,1.0/3.0) + Rj);

    double alpha_param, beta_param;
    if(Zp == 0){
        alpha_param = 0.76+1.93*std::pow(Af,-1.0/3.0);
        beta_param = (1.66*std::pow(Af,-2./3.) - 0.05)/alpha_param;
        }
    else{
        double C_param = 0.0;
        if (Zf >= 50)C_param=-0.10/Ap;
        else if (Zf > 20) {
                C_param =(0.123482-0.00534691*Zf-0.0000610624*Zf*Zf+5.93719*1e-7*Zf*Zf*Zf+1.95687*1e-8*Zf*Zf*Zf*Zf)/Ap;
        }
        /*
        if(Zf<=20)C_param=0.0;
        else if(Zf<=30){int dZ = Zf-20; C_param=0.0 + dZ*0.1*(-0.06);}
        else if(Zf<=40){int dZ = Zf-30; C_param=-0.06 + dZ*0.1*(-0.1+0.06);}
        else if(Zf<=50){C_param = -0.1;}
        else C_param=-0.1;
*/
        alpha_param = 1.0 + C_param;
        beta_param = -Ccoul;
    }

    auto I0 = [](double t) {
        return std::exp(t) - 1.0;
    };

    auto I1 = [](double t, double tx) {
        return (t - tx + 1.0)*std::exp(tx) - t - 1.0;
    };

    auto I2 = [](double s0, double sx) {
        double S = 1.0/std::sqrt(s0);
        double Sx = 1.0/std::sqrt(sx);
        double p1 = S*S*S*( 1.0 + S*S*( 1.5 + 3.75*S*S) );
        double p2 = Sx*Sx*Sx*( 1.0 + Sx*Sx*( 1.5 + 3.75*Sx*Sx) )*std::exp(sx-s0);
        return p1-p2;
    };

    auto I3 = [](double s0, double sx) {
        double s2 = s0*s0;
        double sx2 = sx*sx;
        double S = 1.0/std::sqrt(s0);
        double S2 = S*S;
        double Sx = 1.0/std::sqrt(sx);
        double Sx2 = Sx*Sx;

        double p1 = S *(2.0 + S2 *( 4.0 + S2 *( 13.5 + S2 *( 60.0 + S2 * 325.125 ))));
        double p2 = Sx*Sx2 *(
			 (s2-sx2) + Sx2 *(
					  (1.5*s2+0.5*sx2) + Sx2 *(
								   (3.75*s2+0.25*sx2) + Sx2 *(
											      (12.875*s2+0.625*sx2) + Sx2 *(
															    (59.0625*s2+0.9375*sx2) + Sx2 *(324.8*s2+3.28*sx2))))));

        p2 *= std::exp(sx-s0);
        return p1-p2;
    };

    prefragment proj(Am, Zm);
    prefragment daughter(Af, Zf);

    const double Mp = (Zp==0)?neutron_mass:ame16::get_nuclear_mass(Ap,Zp) * atomic_mass_unit;
    double C = Mp/(nurex::PI*nurex::PI*nurex::hc_transp*nurex::hc_transp);
    double G = (Ap%2 == 0)?1.0:2.0;
    if(Ap==2 && Zp==1)G=3;

    int P = 0;
    if(Zf%2 == 0)P++;
    if((Af-Zf)%2 == 0)P++;
    double P0 = pairing_energy(Af,Zf);
    double Ux = (2.5 + 150.0/Af);
    double Ecrit = Ux +P0;
    double u =  0.05*(Ex - P0);
    double uu = (1-exp(-u))/u;
    daughter.atilda = Af*(0.142  + 0.00917*frdm::get_EFLmic(Af,Zf));
    double a = (daughter.atilda*uu) + ((0.1375 - 8.36E-5*Af)*Af*(1-uu));
    //a = daughter.atilda;
    double T  = 1.0/(std::sqrt(a/Ux) - 1.5/Ux);
    double E0 = Ecrit - T*(std::log(T) - std::log(a)/4.0
	            - 1.25*std::log(Ux) + 2.0*std::sqrt(a*Ux));
    //P = 0;
    //if(Zf%2 == 0)P++;
    //if((Af-Zf)%2 == 0)P++;
    //P0 = 0.5*P*Pa(Am,Zm);
    //double Ecrit_CN = (2.5 + 150.0/Am) + P0;

    double eps_max = Ex - (SB - B) - Ccoul;
    double t = eps_max/T;

    if ( eps_max < Ecrit ) {
        res = (I1(t,t)*T + (beta_param+Ccoul)*I0(t))/std::exp(E0/T);
    }
    else {
        double expE0T = exp(E0/T);
        double sqrt2 = std::sqrt(2.0);

        double tx = Ecrit/T;
        double s0 = 2.0*std::sqrt(a*(eps_max-P0));
        double sx = 2.0*std::sqrt(a*(Ecrit-P0));
        if(s0 > 350.) { s0 = 350.; }
        res = I1(t,tx)*T/expE0T + I3(s0,sx)*exp(s0)/(sqrt2*a) + (beta_param + Ccoul)*(2*sqrt2*I2(s0,sx)*exp(s0));
        //res = I1(t,tx)*T/expE0T + I3(s0,sx)*exp(s0)/(sqrt2*a);
        //printf("2: %lf %lf\n",I1(t,tx), I3(s0,sx));
    }
    double Rb = 0.0;
    if(Ap > 4)
    {
      double Ad = pow(Af,1.0/3.0);
      double Aj = pow(Ap,1.0/3.0);
      Rb = 1.12*(Aj + Ad) - 0.86*((Aj+Ad)/(Aj*Ad))+2.85;
    }
  else if (Ap>1)
    {
      double Ad = pow(Af,1.0/3.0);
      double Aj = pow(Ap, 1.0/3.0);
      Rb=1.5*(Aj+Ad);
    }
  else
    {
      double Ad = pow(Af, 1.0/3.0);
      Rb = 1.5*Ad;
    }
  double GeometricalXS = nurex::PI*Rb*Rb;
  res *=nurex::PI * C*G*GeometricalXS*alpha_param;

  return res;
}


double width_e(prefragment &frag, int Ap, int Zp, double B, double SB, double Ex, double j){
    if(SB>=Ex)return 0.0;
    double res = 0;
    const int Af = frag.A - Ap;
    const int Zf = frag.Z - Zp;
    if(Ap>=Af)return 0.0;
    prefragment daughter(Af,Zf, frag.config);
    const double Rgeom = 1.16 * (pow(Ap,1.0/3.0)+pow(Af,1.0/3.0));
    double Mf = ame16::get_nuclear_mass(Af,Zf) * atomic_mass_unit;
    if(Mf<=0)Mf = Af * atomic_mass_unit;
    const double Mp = (Zp==0)?neutron_mass:ame16::get_nuclear_mass(Ap,Zp) * atomic_mass_unit;
    const double mu = Mf * Mp / (Mf + Mp);
    const double S = SB - B;
    assert(frag.A>Ap);
    assert(Mf>900. && Mp>900. && Mf>Mp);
    const double Ea = Ex - S;
    auto f = [&](double epsilon){
        double Ef = Ea - epsilon;
        double Ecm = epsilon * (Mf-Mp)/Mf;
        double Ra = hc_transp * std::sqrt(0.5/(mu*Ecm));
        assert(std::isfinite(Ra) && (!std::isnan(Ra)));
        //double edif = epsilon - B;
        //double edif = epsilon;
        double rho;
        if(daughter.config.density==level_density_type::GC_GEM){
            rho = level_density_gem(daughter, Ef);
        }
        else if(daughter.config.density==level_density_type::GC_RIPL){
            rho = level_density_ripl(daughter, Ef);
        }
        else if(daughter.config.density==level_density_type::GC_KTUY05){
            rho = level_density_kawano(daughter, Ef);
        }
        else{
            rho = level_density(daughter, Ef, j).first;
        }
        assert( (!std::isnan(rho)) && std::isfinite(rho));
        return (Ra+Rgeom)*(Ra+Rgeom)*epsilon*rho;
    };
    double C = Mp/(nurex::PI*nurex::PI*nurex::hc_transp*nurex::hc_transp);
    double G = (Ap%2 == 0)?1.0:2.0;
    if(Ap==2 && Zp==1)G=3;
    res = integrator_adaptive.integrate(f,B,B+(Ea-B)*0.2,0,1e-3,1);
    res+= integrator_adaptive.integrate(f,B+(Ea-B)*0.2,Ea,0,1e-3,2);
    //res = integrator_adaptive.integrate(f,B,Ea,0,1e-3,1);
    res *= C*G*PI;
    assert( (!std::isnan(res)) && std::isfinite(res));
    return res;
}

double penetration_coefficient(int A, int Ap, double Temp){
    double mu = static_cast<double>((A - Ap) * Ap)/ A;
    const double omega = 4.5;
//    double homega = 197.3287 * omega;
    if(Temp<=0)return 0;
    else{
        return std::pow(10.0,4.e-4*std::pow(Temp/(omega*omega*std::pow(mu,0.25)),-4.3/2.3026));
    }
};

int closest_shell_difference(int n){
    if(n<5){
        return n-2;
    }
    else if(n<11){
        return n-8;
    }
    else if(n<17){
        return n-14;
    }
    else if(n<24){
        return n-20;
    }
    else if(n<39){
        return n-28;
    }
    else if(n<66){
        return n-50;
    }
    else if(n<104){
        return n-82;
    }
    else{
        return n-126;
    }
};

} // end of namespace nurex
