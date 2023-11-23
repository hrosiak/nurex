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
#include "numerics.h"

namespace nurex {

//integrator_1d_type integrator;
GaussLegendreIntegration<integration_1d_order> integratorGL;
GaussLegendreIntegration<8> integratorGL8;
GaussHermiteIntegration<10> integratorGH;
GaussHermiteIntegration2D<6> integratorGH2D;
integrator_adaptive_type integrator_adaptive;

integrator_2d_type integrator2D;

/// interpolation
#ifdef GSL_INTERPOLATION
InterpolatorGSL::InterpolatorGSL(const double *x, const double *y, int num,interpolation_t type){
    acc = gsl_interp_accel_alloc ();

    if(type==cspline)
	spline = gsl_spline_alloc (gsl_interp_cspline, num);
    else
	spline = gsl_spline_alloc (gsl_interp_linear, num);

    gsl_spline_init (spline, x, y, num);
    min= x[0];
    max= x[num-1];

}
InterpolatorGSL::InterpolatorGSL(const std::vector<double>& x, const std::vector<double>& y,interpolation_t type){
    //Interpolator(x.data(),y.data(),x.size());
    acc = gsl_interp_accel_alloc ();
    if(type==cspline)
	spline = gsl_spline_alloc (gsl_interp_cspline, x.size());
    else
	spline = gsl_spline_alloc (gsl_interp_linear, x.size());

    gsl_spline_init (spline, x.data(), y.data(), x.size());
    min= x[0];
    max= x[x.size()-1];
}

InterpolatorGSL::~InterpolatorGSL(){
    gsl_interp_accel_free (acc);
    gsl_spline_free (spline);
}

double InterpolatorGSL::eval(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval(spline, x, acc);
}

double InterpolatorGSL::derivative(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval_deriv (spline, x, acc);
}
#endif

double InterpolatorSpline::eval(double x)const{
    if(x<min)x=min;
    if(x>max)x=max;
    return s(x);
}


double Combination(unsigned int n, unsigned int k){
  double p = 1.0;
  for(unsigned int i=n;i>k;i--){
    p*=i;
  }
  return p/factorial(n-k);
}


} // namespace nurex
