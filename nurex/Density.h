/*
 *  Copyright(C) 2017, Andrej Prochazka, M. Takechi
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

#ifndef DENSITYBASE_H
#define DENSITYBASE_H

#include <memory>
#include <utility>
#include <type_traits>
#include <string>
#include "nurex/numerics.h"
namespace nurex{

/**
 * The enum for available nuclear densities
 */
enum class density_type{fermi, ho, gaussian, dirac, zero, table};

/**
 * Rrms radii of the distribution
 * @tparam T - Density function class
 * @param density 
 * @return @return Rrms radii
 */
template <typename T>
double Rrms(const T& density) noexcept{
    using namespace nurex;
    double res;
    if(density.Norm()==0.0)return 0.0;
    auto f = [&](double x){return x*x*x*x*density.Density(x);};
    double C = 4*PI/density.Norm();
    res = C*integrator_adaptive.integrate(f,0,MAX_INTEGRATION_RADIUS,rrms_integration_precision/C,0,4);
    return sqrt(res);
}

/**
 * helper function to get nuclear density
 * \tparam T - density type
 * \return - pointer to the density type
 Example usage:
 ~~~ {.cpp}
 auto density1 = make_density<DensityFermi>(4.09713,0.555);
 ~~~
 */
template<typename T, typename... Args>
T make_density(Args&&... args){
	if constexpr(std::is_pointer<T>::value){
		return new typename std::remove_pointer<T>::type (std::forward<Args>(args)...);
	}
	else return T(std::forward<Args>(args)...);
}


/**
 * helper function to normalize the nuclear density
 * \tparam T - Density function, ie DensityFermi
 * \param A - mass number for normalization
 */

template<typename T>
double normalize(T& density, double A=0.0){
    double res;
    if(A==0.0){
        A = density.norm;
    }
    auto f = [&](double x){return x*x*density.Density(x);};
    res = 4*PI*integrator_adaptive.integrate(f,0,MAX_INTEGRATION_RADIUS,rrms_integration_precision/(4*PI),0,4);
    density.rho0=A*density.rho0/res;
    density.norm = A;
    return A;
    }

//************* DensityFermi ***************//
/** 
 *  Fermi function
 */
class DensityFermi{
	public:
	DensityFermi(double radius, double diffuseness, double w=0.0, double normalization=1.0);
	double Density(double r)const noexcept;
    double Normalize(double a=0){return normalize(*this,a);}
	double GetParameter(int i) const noexcept;
    double Norm()const noexcept {return norm;}

    void SetParameters(double radius, double diffuseness, double w=0.0);
    double GetRadiusParameter(){return r0;}
    double GetDiffusenessParameter(){return c;}
	/// returns normalization parameter
    double GetRho0(){return rho0;}

    double operator()(double r)const noexcept{return Density(r);}

    constexpr static int nparameters = 3;
    constexpr static density_type type = density_type::fermi;
	double r0=1;
	double c=0.5;
	double rho0=1;
    double norm=1;
	double w=0.0;
	};


//************* DensityHO ***************//
/** 
 *  Harmonic oscillator type density function
 */
class DensityHO{
	public:
	DensityHO(double radius, double _a, double normalization=1);
	double Density(double r)const noexcept;
    double Normalize(double a=0){return normalize(*this,a);}
	double GetParameter(int i) const noexcept;
	void SetParameters(double radius, double width);
    double Norm()const noexcept {return norm;}
	/// returns Radius parameter \f$r_0\f$
    double GetRadiusParameter(){return r0;}
	/// returns Difuseness parameter a
    double GetDiffusenessParameter(){return a;}
	/// returns normalization parameter
    double GetRho0(){return rho0;}
    double operator()(double r)const noexcept{return Density(r);}
    constexpr static int nparameters = 2;
    constexpr static density_type type = density_type::ho;
    double r0=0;
	double a=1;
	double rho0=1;
    double norm=1;
	};

//************* DensityGaussian ***************//
/** 
 *  Gaussian type density function
 */
class DensityGaussian{
public:
DensityGaussian(double width, double normalization=1);
double Density(double r)const noexcept;
double Normalize(double A=0);
double GetParameter(int i) const noexcept;
void SetParameters(double width);
double Norm()const noexcept {return norm;}
double GetWidth(){return width;}
/// returns normalization parameter
double GetRho0(){return rho0;}

double operator()(double r)const noexcept{return Density(r);}
constexpr static int nparameters = 1;
constexpr static density_type type = density_type::gaussian;
double width=1;
double rho0=1;
double norm=1;
};

//************* DensityDirac ***************//
/** 
 *  Dirac density function
 */
class DensityDirac{
	public:
    DensityDirac(double normalization=1):norm(normalization){}
    double Density(double)const noexcept{return 0.0;}
    double Normalize(double A=1){norm=A;return norm;}
	double GetParameter(int) const noexcept {return 0.0;}
	double Rrms(){return 0.0;}
    double Norm()const noexcept {return norm;}
    double operator()(double r)const noexcept{return Density(r);}
    constexpr static int nparameters = 0;
    constexpr static density_type type = density_type::dirac;
    double norm=1;
	};

/** 
 *  Zero Density, usefull for proton or neutron distributions
 */
//************* DensityZero ***************//
class DensityZero{
	public:
    DensityZero(){}
    static constexpr double Density(double )noexcept{return 0.0;}
    static constexpr double Normalize(double){return 0.0;}
    static constexpr double GetParameter(int) noexcept {return 0.0;}
    static constexpr double Norm()noexcept {return 0.0;}
    constexpr double operator()(double ) const noexcept{return 0.0;}
    constexpr static int nparameters = 0;
    constexpr static density_type type = density_type::zero;
    };

/** 
 *  tabulated density
 */
class DensityTable{
    public:
    DensityTable(std::vector<double> b, std::vector<double> d, double normalization = 1.0);
    ~DensityTable()
    ;
    double Density(double r)const noexcept;
	double Normalize(double A=0);
    double GetParameter(int) const noexcept {return 0.0;}
    double Norm()const noexcept {return norm;}
    double operator()(double r)const noexcept{return Density(r);}
	/// returns normalization parameter
    double GetRho0(){return rho0;}
    constexpr static int nparameters = 0;
    constexpr static density_type type = density_type::table;
    double norm=1.0;
    double rho0=1.0;
    private:
    Interpolator density;
    double max=0.0;
};

/**
 * DensityType - conceptual class which can hold any Density class with required functions
 * It can store class with following functions and variables:
 *   * double Density(double)
 *   * double Normalize(double)
 *   * double GetParameter(int i) const
 *   * double Norm()const noexcept
 *   * int nparameters  // number of used parameters
 *   * double Norm()const noexcept;
 *   * density_typpe type  // density type identifier
 */
class DensityType{
    struct concept_t{
        virtual ~concept_t(){}
        virtual double Density(double) const noexcept=0;
        virtual double Normalize(double a=0)=0;
        virtual double GetParameter(int i) const noexcept=0;
        virtual int NumberOfParameters() const noexcept=0;
        virtual double Norm()const noexcept = 0;
        virtual density_type type() const noexcept=0;
        virtual std::unique_ptr<concept_t> Clone() const=0;
    };

    template<typename T>
    struct model_t final:public concept_t{
        model_t() = default;
        model_t(const T& v) : m_data(v) {}
        model_t(T&& v) : m_data(std::move(v)) {}
        double Density(double r)const noexcept override{return m_data.Density(r);}
        double Normalize(double a=0)override{return m_data.Normalize(a);}
        double GetParameter(int i)const noexcept override{return m_data.GetParameter(i);}
        int NumberOfParameters() const noexcept override {return m_data.nparameters;}
        double Norm() const noexcept override{return m_data.Norm();}
        density_type type() const noexcept override {return m_data.type;}
        std::unique_ptr<concept_t> Clone() const override {return std::make_unique<model_t>(*this);}
        T m_data;

        static_assert(!std::is_const_v<T>, "Shouldn't be const here");
        static_assert(!std::is_pointer_v<T>, "Shouldn't be pointer here");
    };
public:
    DensityType()=default;
    DensityType(const DensityType& other):object(other.Clone()){}
    DensityType(DensityType&&) = default;

    template<typename T>
    DensityType(T x):object(std::make_unique<model_t<T>>(std::move(x))){}
    DensityType& operator=(const DensityType& other){
        object = other.Clone();
        return *this;
    }
    DensityType& operator=(DensityType&&) = default;
    
    /// root-mean-squared radius 
    double Rrms()const {if(!object)return -1.0;if(rrms>=0.0)return rrms; rrms = nurex::Rrms(*object); return rrms;}

    /// return  density at radius r
    double Density(double r)const {return object->Density(r);}

    /// Normalize density to a 
    double Normalize(double a=0){return object->Normalize(a);}

    /// get i-th parameter
    double GetParameter(int i)const {return object->GetParameter(i);}

    /// return number of parameters
    int NumberOfParameters() const {return object->NumberOfParameters();}

    /// return normalization
    double Norm() const {return object->Norm();}

    /// return density type identifier
    density_type type() const {return object->type();}

    /// get pointer to the copy of this DensityType  
    std::unique_ptr<concept_t> Clone() const {return object->Clone();}

    /// check if DensityType contains some Density
    explicit operator bool()const{return (object)?true:false;}

private:
    std::unique_ptr<concept_t> object;
    mutable double rrms=-1.0;

};


#ifndef NO_FILESYSTEM
DensityTable density_from_file(const char *filename, double normalization=1.0);
#endif

/**
 * compares two density types
 */
inline bool operator==(const DensityType &d1, const DensityType &d2){
    if( (d1.type() != d2.type()) || (d1.Norm()!=d2.Norm())){
        return false;
    }
    else{
        for(int i=0;i<d1.NumberOfParameters();i++){
            if(d1.GetParameter(i)!=d2.GetParameter(i)){
                return false;
            }
        }
    }
    return true;
}

inline bool operator!=(const DensityType &d1, const DensityType &d2){
    return !(d1==d2);
}

/**
 * converts density type to string
 */
std::string density_type_to_string(const density_type dt);

/**
 * converts DensityType class to string
 */
inline std::string density_type_to_string(const DensityType& d){
	return density_type_to_string(d.type());
}

} // end of namespace
#endif
