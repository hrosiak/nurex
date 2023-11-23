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

#include "ModelUtils.h"

namespace nurex {


double z_integral(const DensityType& density, double b)
{
    double res;
    double rmax = max_rrms_factor * density.Rrms();

    auto f = [&](double _z) {
        return density.Density(sqrt((b * b) + (_z * _z)));
    };

    res = 2 * integratorGL(f, 0, rmax);
    //res = 2 * integrator_adaptive.integrate(f, 0, rmax,0,zintegration_precision,2);
    return res;
}

Functional ZIntegrate(const DensityType& nd)
{
    // check for special cases
    if (nd.Rrms() == 0.0 && nd.Norm() == 0.0) {
        return ConstantFunction(0.0);
    }

    if (nd.type() == density_type::dirac) {
        return DiracFunction(nd.Norm());
    }

    double rmax = max_rrms_factor * nd.Rrms();
    LinearVArray<z_integration_b_steps> x(0, rmax);
    std::array<double,z_integration_b_steps> zi;

    for(unsigned int i = 0; i<z_integration_b_steps; i++){
        zi[i] = z_integral(nd, x[i]);
    }
    // this is to add last points as 0.0
    zi[z_integration_b_steps - 1] = 0.0;
    return InterpolatorSplineT<LinearVArray<z_integration_b_steps>>(x, zi);
}

Functional ZIntegrateRange(const DensityType& nd, double beta)
{
    // check for special cases
    if (nd.Rrms() == 0.0 && nd.Norm() == 0.0) {
        return ConstantFunction(0.0);
    }

    if (nd.type() == density_type::dirac) {
        return DiracFunction(nd.Norm());
    }

    double rmax = max_rrms_factor * nd.Rrms();
    const int n = z_range_integration_b_steps;
    LinearVArray<n> x(0, rmax);
    std::array<double,n> zi;
    //double radius;
    auto f_profile_rho = [&](double x, double y) {
            //double r2 = distance2(0, radius, x, y);
            //return z_integral(nd,distance(x, y, 0.0, 0.0)) * finite_range2(r2, beta);
            return z_integral(nd,distance(x, y, 0.0, 0.0));
        };
    for(unsigned int i = 0; i<n; i++){
        zi[i] = 0.5*integratorGH2D(f_profile_rho,0, beta, x[i] ,beta)/(PI*beta*beta);
        //radius = x[i];
        //zi[i] = integrator2D(f_profile_rho,-max_finite_range_factor*beta,max_finite_range_factor*beta, x[i] - max_finite_range_factor*beta, x[i] + max_finite_range_factor*beta);
    }
    return InterpolatorSplineT<LinearVArray<n>>(x, zi);
}

double z_fermi_energy(const DensityType& d, double mass, double b){
    double rmax = max_rrms_factor * d.Rrms();

    auto f = [&](double z) {
        double r = sqrt((b*b)+(z*z));
        return d.Density(r)*fermi_energy(d.Density(r),mass);
    };
    auto norm = [&](double z) {
        double r = sqrt((b*b)+(z*z));
        return d.Density(r);
    };

    return integratorGL(f,0,rmax)/integratorGL(norm,0,rmax);
}

Functional ZIntegrate_Fermi_Energy(const DensityType& nd, double mass)
{
    // check for special cases
    if (nd.Rrms() == 0.0 || nd.Norm() == 0.0) {
        return ConstantFunction(0.0);
    }
    double rmax = max_rrms_factor * nd.Rrms();
    LinearVArray<z_integration_b_steps> x(0, rmax);
    std::array<double,z_integration_b_steps> zi;

    for(unsigned int i = 0; i<z_integration_b_steps; i++){
        zi[i] = z_fermi_energy(nd, mass, x[i]);
    }
    // this is to add last points as 0.0
    zi[z_integration_b_steps - 1] = 0.0;

    return InterpolatorSplineT<LinearVArray<z_integration_b_steps>>(x, zi);
}

double z_fermi_momentum(const DensityType& d, double b){
    double rmax = max_rrms_factor * d.Rrms();

    auto f = [&](double z) {
        double r = sqrt((b*b)+(z*z));
        return d.Density(r)*fermi_momentum(d.Density(r));
    };
    auto norm = [&](double z) {
        double r = sqrt((b*b)+(z*z));
        return d.Density(r);
    };

    return integratorGL(f,0,rmax)/integratorGL(norm,0,rmax);
}

Functional ZIntegrate_Fermi_Momentum(const DensityType& nd)
{
    // check for special cases
    if (nd.Rrms() == 0.0 || nd.Norm() == 0.0) {
        return ConstantFunction(0.0);
    }
    double rmax = max_rrms_factor * nd.Rrms();
    LinearVArray<z_range_integration_b_steps> x(0, rmax);
    std::array<double,z_range_integration_b_steps> zi;

    for(unsigned int i = 0; i<z_range_integration_b_steps; i++){
        zi[i] = z_fermi_momentum(nd, x[i]);
    }
    // this is to add last points as 0.0
    zi[z_range_integration_b_steps - 1] = 0.0;

    return InterpolatorSplineT<LinearVArray<z_range_integration_b_steps>>(x, zi);
}


}