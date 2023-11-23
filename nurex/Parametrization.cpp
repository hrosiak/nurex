#include "nurex/Parametrization.h"
#include "nurex/Nucleus.h"
#include "nurex/Config.h"
#include "nurex/Utils.h"
#include "cmath"
#include <iostream>
using std::pow;

namespace nurex{

double SigmaR_Kox(int Ap, int Zp, double E, int At, int Zt){
    constexpr double rc = 1.3;
    constexpr double r0 = 1.1;
    constexpr double a = 1.85;
    constexpr double c1 = 2.0-(10.0/(1.5*1.5*1.5*1.5*1.5));
    double Ap13 = pow(Ap,1.0/3.0);
    double At13 = pow(At,1.0/3.0);
    double D = 5.0*(At-2*Zt)*Zp/(Ap*At);
    double Bc = Zp*Zt/(rc*(Ap13+At13));
    double logE = std::log10(E);
    double c = 0;
    if(logE < 1.5){
        c = c1*std::pow(logE/1.5,3);
    }
    else{
        c = (-10.0/std::pow(logE,5)) + 2;
    }
    double Rs = r0 * ((a*Ap13*At13)/(Ap13+At13)-c)+D;
    double Rv = r0 * (Ap13 + At13);
    double Ri = Rv + Rs;
    //double Ecm = Ecm_from_T_relativistic(E,Ap,At);
    double Ecm = Ecm_from_T(E,Ap,At);
    return 10.0*PI*Ri*Ri*(1-(Bc/Ecm));
    
    }
double SigmaR_Kox(const Nucleus &projectile, double E, const Nucleus &target){
    return SigmaR_Kox(projectile.A(), projectile.Z(), E, target.A(), target.Z());
    }
} //end of nurex namespace
