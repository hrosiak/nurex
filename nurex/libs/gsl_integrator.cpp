#include "gsl_integrator.h"
#include "math.h"
#include "gsl/gsl_errno.h"

namespace integrators{

double funcwrapper3(double x, void *_c){
    std::function<double(double)> *f = (std::function<double(double)> *)_c;
    return (*f)(x);
    }

IntegratorGSL::IntegratorGSL(bool adapt):adaptive(adapt){
        gsl_set_error_handler_off();
        if(adaptive){
            w=gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
        }
    };

IntegratorGSL::~IntegratorGSL(){
        if(adaptive){
            gsl_integration_workspace_free(w);
        }
    };

double IntegratorGSL::integrate(std::function<double(double)> f, double _min, double  _max, double absprecision, double prec){
        gsl_function F;
        F.function = funcwrapper3;
        F.params = &f;

        min = _min;
        max = _max;
        size_t num;
        if(adaptive){
            #ifdef USE_THREADS
            std::lock_guard<std::mutex> lock(integration_mutex);
            #endif
            gsl_integration_qag(&F,_min,_max,absprecision,prec,GSL_WORKSPACE_SIZE,6,w,&result,&error);
        }
        else
            gsl_integration_qng(&F,_min,_max,absprecision,prec,&result,&error,&num);
        return result;
    };


IntegratorGSL2::IntegratorGSL2(std::function<double(double)> *f, bool adapt):_f(f),adaptive(adapt){
        gsl_set_error_handler_off();
        if(adaptive){
            w=gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
        }
    };

IntegratorGSL2::~IntegratorGSL2(){
        if(adaptive){
            gsl_integration_workspace_free(w);
        }
    };

double IntegratorGSL2D::operator()(std::function<double(double,double)> f, double min1, double max1, double min2, double max2){\
    double par;
    double res;

    auto f1 = [&](double y){
            return f(par,y);
        };

    auto f2 = [&](double x){
        par = x;
        return i1->integrate(f1,min2,max2,1e-2);
    };
    res = i2->integrate(f2, min1,max1,1e-2);
    return res;
};

double IntegratorGSL2::integrate(double _min, double  _max, double prec){
        gsl_function F;
        F.function = funcwrapper3;
        //F.function = static_cast<gsl_function *>(&fw);
        F.params = _f;

        min = _min;
        max = _max;
        size_t num;
        if(adaptive){
            #ifdef USE_THREADS
            std::lock_guard<std::mutex> lock(integration_mutex);
            #endif
            gsl_integration_qag(&F,_min,_max,1e-6,prec,GSL_WORKSPACE_SIZE,6,w,&result,&error);
        }
        else
            gsl_integration_qng(&F,_min,_max,1e-6,prec,&result,&error,&num);
        return result;
    };


double Integrate2D(std::function<double(double,double)> f, double min1, double max1, double min2, double max2){
    double par;
    double res;

    auto f1 = [&](double y){
            return f(par,y);
        };
    std::function<double(double)> _f1 = f1;
    IntegratorGSL2 i1(&(_f1),true);

    auto f2 = [&](double x){
        par = x;
        return i1.integrate(min2,max2,1e-2);
    };

    std::function<double(double)> _f2 = f2;
    IntegratorGSL2 i2(&(_f2),true);
    res = i2.integrate(min1,max1,1e-2);

    return res;
}

double Integrate2Dr(std::function<double(double,double)> &f, double min1, double max1, double min2, double max2){
    double par;

    auto f1 = [&](double y){
            return f(par,y);
        };
    std::function<double(double)> _f1 = f1;
    IntegratorGSL2 i1(&(_f1),true);

    auto f2 = [&](double x){
        par = x;
        return i1.integrate(min2,max2,1e-2);
    };

    std::function<double(double)> _f2 = f2;
    IntegratorGSL2 i2(&(_f2),true);
    return i2.integrate(min1,max1,1e-2);
}

} //end of sigmar namespace
