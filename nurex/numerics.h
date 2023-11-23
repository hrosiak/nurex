/*
 *  Copyright(C) 2017, Andrej Prochazka
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <cmath>
#include <array>
#include <functional>
#include <memory>
#include <type_traits>
#include "nurex/Config.h"
#include "nurex/libs/glq_integrator.h"
#include "nurex/libs/gkq_integrator.h"
#include "nurex/libs/ghq_integrator.h"
#include "nurex/libs/spline.h"

#ifdef GSL_INTEGRATION
#include "nurex/gsl_integrator.h"
#endif

#ifdef USE_THREADS
#include <mutex>
#endif

using integrators::GaussLegendreIntegration;
using integrators::GaussLegendreIntegration2D;
using integrators::GaussLegendreIntegration2DA;
using integrators::GaussHermiteIntegration;
using integrators::GaussHermiteIntegration2D;
using integrators::GaussKronrodIntegration;

#ifdef GSL_INTEGRATION
using integrators::IntegratorGSL;
using integrators::IntegratorGSL2D;
#endif

namespace nurex{
#ifdef GSL_INTEGRATION
using  integrator_2d_type = IntegratorGSL2D;
using  integrator_fm_type = IntegratorGSL;
using  integrator_adaptive_type =  IntegratorGSL;
#else
using  integrator_2d_type = GaussLegendreIntegration2DA<integration_2d_order,8>;
using  integrator_fm_type = GaussKronrodIntegration<21>;
using  integrator_adaptive_type =  GaussKronrodIntegration<21>;
#endif

extern integrator_adaptive_type integrator_adaptive;
extern GaussLegendreIntegration<integration_1d_order> integratorGL;
extern GaussLegendreIntegration<8> integratorGL8;
extern GaussHermiteIntegration<10> integratorGH;
extern GaussHermiteIntegration2D<6> integratorGH2D;
extern integrator_2d_type integrator2D;


/// interpolation ////

#define IS_DIRAC(x,y) is_ftype<DiracFunction>((x)) || is_ftype<DiracFunction>((y))

template<class T, class=void>
struct has_minmax : std::false_type{};
template<class T>
struct has_minmax<T, std::void_t<decltype(std::declval<T>().get_min()),decltype(std::declval<T>().get_max())>> : std::true_type{};

template<class T, class=void>
struct has_eval : std::false_type{};
template<class T>
struct has_eval<T, std::void_t<decltype(std::declval<T>().eval(0.))>> : std::true_type{};


/**
 * class to store 1-d function like objects
 */ 
class Functional{
    struct concept_t{
        virtual ~concept_t(){}
        virtual double eval(double) const = 0;
        virtual double min()const = 0;
        virtual double max()const = 0;
    };

    template<typename T>
    struct model_t final:public concept_t{
        static_assert(!std::is_pointer_v<T>, "Shouldn't be pointer here");
        model_t() = default;
        model_t(const T& v) : m_data(v) {}
        model_t(T&& v) : m_data(std::move(v)) {}
        double eval(double x) const override final {
            if constexpr (has_eval<T>::value) return m_data.eval(x);
            else return m_data(x);
            }
        double min()const override final {
            if constexpr(has_minmax<T>::value)return m_data.get_min();
            else return 0.0;}
        double max()const override final {
            if constexpr(has_minmax<T>::value)return m_data.get_max();
            else return 0.0;}
        T m_data;
    };

public:
    Functional()=default;
    Functional(const Functional&)=delete;
    Functional(Functional&&) = default;
    Functional& operator=(Functional&&) = default;
    Functional& operator=(const Functional&) = delete;

    template<typename T>
    Functional(T x):object(std::make_unique<model_t<T>>(std::move(x))){}

    template<typename T>
    bool is_type() const {return dynamic_cast<model_t<T>*>( object.get() ) != nullptr;}


    double eval(double x) const {return object->eval(x);}
    double operator()(double r) const {return eval(r);}
    double min() const {return object->min();}
    double max() const {return object->max();}
    private:
        std::unique_ptr<concept_t> object;
};


/// The enum for interpolation types
enum interpolation_t {cspline, linear};

#ifdef GSL_INTERPOLATION
/// Interpolation class, to store interpolated values
class InterpolatorGSL{
        public:
        InterpolatorGSL(const double *x, const double *y, int num,interpolation_t type=cspline);
        InterpolatorGSL(const std::vector<double>& x, const std::vector<double>& y,interpolation_t type=cspline);
        ~InterpolatorGSL();
        double operator()(double x){return eval(x);};
        double eval(double x);
        double derivative(double x);
        double get_min(){return min;};
		double get_max(){return max;};
		double integral(double lo, double hi);

        private:
        double min=0;
        double max=0;
        gsl_interp_accel *acc;
        gsl_spline *spline;
	};
#endif

class InterpolatorSpline{
    public:
    InterpolatorSpline(const std::vector<double>& x, const std::vector<double>& y,bool cspline = true):s(cspline_vector(x,y,cspline)){
        //s = std::move();
        min=x[0];max=x[x.size()-1];
    }
        double operator()(double x)const{return eval(x);}
        double eval(double x) const;
        double get_min()const{return min;}
        double get_max()const{return max;}
    private:
        cspline_vector s;
        double min=0;
        double max=0;
};

template<typename T>
class InterpolatorSplineT{
    public:
        InterpolatorSplineT(T x, const std::array<double,T::size()> y):s(std::move(cspline_special<T>(x,y))){
            min = x[0];
            max = x[x.size()-1];
        }
        double operator()(double x){return eval(x);}
        double eval(double x) const {
            if(x<min)x=min;
            if(x>max)x=max;
            return s.evaluate(x);
        }
        double get_min()const{return min;}
        double get_max()const{return max;}
    private:
        cspline_special<T> s;
        double min=0;
        double max=0;
};


#ifdef GSL_INTERPOLATION
using Interpolator = InterpolatorGSL;
#else
using Interpolator = InterpolatorSpline;
#endif

class ConstantFunction{
public:
	ConstantFunction(double v=0.0):val(v){}
    double eval(double)const{return val;}
    double derivative(double)const{return 0.0;}
private:
	double val = 0.0;
};

class DiracFunction{
public:
	DiracFunction(double v=1.0):val(v){}
    double eval(double x)const{return (x==0.0)?val:0.0;}
    double derivative(double)const{return 0.0;}
private:
	double val = 1.0;
};

///// helper functions //////
// return vector with lineary spaced elements from a to b, num is number of elements
inline std::vector<double> logspace_vector(double a, double b, unsigned int num){
    std::vector<double> res;
    if(num>=2 && a<b && a>0.0){
        double emin = std::log10(a);
        double emax = std::log10(b);
	    res.resize(num);
	    double step = (emax-emin)/(num-1);
	    for(unsigned int i=0;i<(num-1);i++){
    	    res[i]=std::pow(10,emin+(i*step));
	        }
	    res[num-1] = std::pow(10,emax);
	    }
    return res;
    }

/**
 * return vector with lineary spaced elements from a to b, num is number of elements
 */ 
inline std::vector<double> linspace_vector(double a, double b, unsigned int num){
    std::vector<double> res;
    if(num>=2 && a<b){
	res.resize(num);
	double step = (b-a)/(num-1);
	for(unsigned int i=0;i<(num-1);i++){
	    res[i]=a+(i*step);
	    }
	res[num-1] = b;
	}
    return res;
    }

// return vector with lineary spaced elements from a to b, num is number of elements
template <int N>
std::array<double,N> linspace_array(double a, double b){
    static_assert (N>2, "N must be more than 2");
    std::array<double,N> res;
    if(a<b){
    double step = (b-a)/N;
    for(unsigned int i=0;i<(N-1);i++){
        res[i]=a+(i*step);
        }
    res[N-1] = b;
    }
    return res;
    }

// factorial
inline constexpr double factorial(unsigned int n){
    assert(n>=0 && n<100);    
    return (n<=1)?1 : (n*factorial(n-1));
};

template<typename F>
double crossing(F& f, double a, double b, double off = 0.0,double eps = 1e-3){
    if(f(a)<off){
        return a-1.0;
    }
    double c=a;
    while((b-a)>=eps){
        c=(a+b)*0.5;
        if(f(c)>off){
            a = c;
            }
        if(f(c)<=off){
            b = c;
        }
    }
    return c;
}

template<typename F>
double bisection(F& f, double a, double b, double eps = 1e-3){
    if(f(a)*f(b)>=0){
        return a-1.0;
    }
    double c=a;
    while((b-a)>=eps){
        c=(a+b)*0.5;
        if(f(c)==0.0){
            break;
            }
        else if(f(c)*f(a)<0){
            b = c;
        }
        else{
            a = c;
        }
    }
    return c;
}

template<typename F>
double bisection_offset(F& f, double a, double b, double base=0.0, double eps = 1e-3){
    if((f(a)-base)*(f(b)-base)>=0){
        return a-1.0;
    }
    double c=a;
    while((b-a)>=eps){
        c=(a+b)*0.5;
        if((f(c)-base)==0.0){
            break;
            }
        else if((f(c)-base)*(f(a)-base)<0){
            b = c;
        }
        else{
            a = c;
        }
    }
    return c;
}

template<typename F>
double golden_search(const F &f, double a, double b, double tol=1e-3){
        const double gr = (sqrt(5) + 1) / 2.0;
    	double c = b - (b - a) / gr;
    	double d = a + (b - a) / gr;

        while (std::abs(c - d) > tol){
        	if(f(c) > f(d)){
              b = d;
            }
        	else{
              a = c;
            }

        	c = b - (b - a) / gr;
        	d = a + (b - a) / gr;
    }
    return (b + a) / 2.0;
}

double Combination(unsigned int n, unsigned int k);

std::vector<std::pair<int,double>> get_discrete_gaussian(double mean, double sigma);

}//end of nurex namespace

#endif

