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

#include "nurex/Density.h"
#include "nurex/Utils.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
namespace nurex{

/// DensityFermi Constructor
/// by default the distribution is normalized to 1.0
/// @param radius - R0 radius parameter in fm
/// @param diffusion - c diffusion parameter in fm
/// @param normalization - initial normalization factor
DensityFermi::DensityFermi(double radius, double diffusion, double _w, double normalization):r0(radius),c(diffusion),w(_w){
	Normalize(normalization);
	}

/// Density(r)
/// returns the density at radius r in fm^-3
double DensityFermi::Density(double r)const noexcept{
    double f = 1.0;
    if(w!=0.0){
      f = 1.0 + (w*r*r/(r0*r0));
    }
	return rho0*f/(1+exp((fabs(r)-r0)/c));
	}

/// SetParameter - set the radius and diffuseness parameters
void DensityFermi::SetParameters(double radius, double diffuseness, double _w){
    r0 = radius;
    c = diffuseness;
    w=_w;
    Normalize(Norm());
}

double DensityFermi::GetParameter(int i) const noexcept
{
	if(i==0){
		return r0;
		}
	else if(i==1){
		return c;
	}
	else if(i==2){
		return w;
	}
	else{
		return 0.0;
		}
}

//////////////// Density HO /////////////////////
/// DensityHO Constructor
/// by default the distribution is normalized to 1.0
/// @param radius - R0 radius parameter in fm
/// @param diffusion - a diffusion parameter in fm
DensityHO::DensityHO(double radius, double _a, double normalization):r0(radius),a(_a){
	Normalize(normalization);
	}

double DensityHO::Density(double r)const noexcept{
	double r2 = r*r/(r0*r0);
	return rho0*(1+a*r2)*exp(-r2);
	}

/// SetParameter - set the radius and width parameters
void DensityHO::SetParameters(double radius, double width){
    r0 = radius;
    a = width;
    Normalize(Norm());
}

double DensityHO::GetParameter(int i) const noexcept
{
	if(i==0){
		return r0;
		}
	else if(i==1){
		return a;
	}
	else{
		return 0.0;
		}
}

//////////////// Density Gaussian /////////////////////
/// DensityGaussian Constructor
/// by default the distribution is normalized to 1.0
/// @param width - Gaussian width parameter in fm
DensityGaussian::DensityGaussian(double width, double normalization):width(width){
	Normalize(normalization);
	}

double DensityGaussian::Density(double r)const noexcept{
	double r2 = r*r/(width*width);
	return rho0*exp(-0.5*r2);
	}


double DensityGaussian::Normalize(double A){
    return normalize(*this, A);
	//rho0 = A/(pow(width,3.0)*pow(2*PI,3.0/2.0));
    //return rho0;
	}

/// SetParameter - set the radius and width parameters
void DensityGaussian::SetParameters(double w){
    width = w;
    Normalize(Norm());
}

double DensityGaussian::GetParameter(int i) const noexcept
{
	if(i==0){
		return width;
		}
	else{
		return 0.0;
		}
}

//////////////// DensityTable /////////////////////
/// Constructor
/// @param b - vector of b values
/// @param d - vector of corresponding densities
/// @param normalization - optional normalization, default is 1.0
DensityTable::DensityTable(std::vector<double> b, std::vector<double> d, double normalization):density(Interpolator(b,d)),max(b.back()){
    normalize(*this,normalization);
  }

/// Destructor
DensityTable::~DensityTable(){
    //delete density;
}

/// Density
double DensityTable::Density(double r)const noexcept{
    return (r>max)?0.0:rho0*density.eval(r);
}

///Normalize
double DensityTable::Normalize(double A){
    return normalize(*this,A);
  }


#ifndef NO_FILESYSTEM
DensityTable density_from_file(const char *filename, double normalization){
    std::string line;
    std::ifstream file;
    file.open(filename);
    if(!file.is_open()){
        throw std::invalid_argument("Failed to open the file");
    }
    double _x,_y;
    std::vector<double> x,y;
    while(file.good() && file>>_x && file>>_y ){
        if(x.size()>0 && _x<=x.back()){
            printf("not increasing radius:  %lf, prev.:%lf \n",_x,x.back());
        }
        if(_y<0){
            printf("density than 0:  %lf, cc: %lf \n",_x,_y);
        }
        x.push_back(_x);
    if(_y<1e-11){
        y.push_back(0.0);
        break;
    }
    y.push_back(_y);
    }
    return DensityTable(x,y, normalization);
  }
#endif


std::string density_type_to_string(const density_type dt){
	std::string s;
	if(dt == density_type::fermi)s="fermi";
	else if(dt == density_type::ho)s="ho";
	else if(dt == density_type::dirac)s="dirac";
	else if(dt == density_type::zero)s="zero";
	else if(dt == density_type::table)s="table";
	else if(dt == density_type::gaussian)s="gaussian";
	else s="none";
	return s;
}

}
