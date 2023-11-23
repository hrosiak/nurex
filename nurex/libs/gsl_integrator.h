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
#ifndef GSL_INTEGRATOR_H
#define GSL_INTEGRATOR_H

#include <functional>
#include <gsl/gsl_spline.h>
#include "gsl/gsl_integration.h"
#include <vector>

#ifdef USE_THREADS
#include <mutex>
#endif

namespace integrators{

#define GSL_WORKSPACE_SIZE 100

/// helper class to integrate functions using the GSL library
class IntegratorGSL{
	public:
	IntegratorGSL(bool adapt=true);
	~IntegratorGSL();

	double integrate(std::function<double(double)> f, double min, double  max, double absprecision = 1e-3, double precision=1e-3);

	private:
	gsl_integration_workspace *w;
	bool adaptive = true;
	double error=0.0;
	double precision=0.0;
	double result=0.0;
	double min=0.0;
	double max=0.0;
	#ifdef USE_THREADS
	std::mutex integration_mutex;
	#endif
    };

class IntegratorGSL2{
	public:
	IntegratorGSL2(std::function<double(double)> *f, bool adaptive=true);
	~IntegratorGSL2();

	double integrate(double min, double  max, double precision=1e-3);
	
	double error;
	double precision;
	double result;

	private:
	std::function<double(double)> *_f;
	gsl_integration_workspace *w;
	bool adaptive = true;
	double min;
	double max;
	#ifdef USE_THREADS
	std::mutex integration_mutex;
	#endif
    };

class IntegratorGSL2D{
	private:
	IntegratorGSL *i1,*i2;
	public:
	IntegratorGSL2D(){i1 = new IntegratorGSL(false);i2 = new IntegratorGSL(false);}
	~IntegratorGSL2D(){delete i1; delete i2;}
	double operator()(std::function<double(double,double)> f, double min1, double max1, double min2, double max2);
};

/// function to integrate 2D function using the GSL library
double Integrate2D(std::function<double(double,double)> f, double min1, double max1, double min2, double max2);
double Integrate2Dr(std::function<double(double,double)> &f, double min1, double max1, double min2, double max2);

}
#endif
