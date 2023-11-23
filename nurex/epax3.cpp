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

#include <math.h>
#include <stdio.h>
namespace epax {

/// epax3 coefficients
const double s[2] = {0.27,1.8};    // scaling in barn
const double p[2] = {-1.731,-0.01399};
const double y[2] = {0.75,0.1};
const double d1 =-1.087;
const double d2 = 3.047E-2;
const double d3 = 2.135E-4;
const double d4 = 71.35;
const double d5 =-25.0;
const double d6 = 0.8;

const double n1 = 0.4;
const double n2 = 0.6;

const double q1 = -10.25;
const double q2 =  10.25;

const double r[7] = {2.78, -0.015, 3.2e-5, 0.0412, 0.124, 30.0, 0.85};

const double Un = 1.65;
const double Up = 2.1;

const double l1 = 1.2;
const double l2 = 0.647;

const double b1 = 2.3e-3;
const double b2= 2.4;

double epax3(int Ap, int Zp, int At, int Zt, int A, int Z){  
  double Zprob,Zbeta, Zbeta_p;

  double R,R0,zp_diff,delta,delta_m;
  //double expo_fac,offset,dzprob;
  double expo,fract,yield_a;
  double norm;
  double result;
  result = -1.0;

  double Ap13 = pow(Ap,1.0/3.0);
  double At13 = pow(At,1.0/3.0);
  
    /* calculate mass yield: */

  /* slope parameter */
  double P = exp(p[0]+p[1]*Ap) ;
  double S = s[0]*P*(Ap13 +At13-s[1]) ;
  yield_a = S * exp(-P * (Ap - A)) ;
  double mass_ratio = (double)A/Ap;

  if ( mass_ratio > y[0]) {
    yield_a *= exp(y[1]*Ap*pow( ((double)A/Ap)-y[0],2)) ;
  }

  // centroid Zprob
  Zbeta = (double)A/(1.98+0.0155*pow(A,2.0/3.0)) ;  // Eq (7)
  Zbeta_p = (double)Ap/(1.98+0.0155*pow(Ap,2.0/3.0));

  if(A > d4) {  // Eq (8)
     delta = d1+d2*A ;
  }
  else
  {
     delta = d3*A*A ; 
  } 

  if(mass_ratio > d6) {  // Eq. (9)
    delta = delta*( 1 + (d5*pow(mass_ratio-d6,2)) );
  }

  zp_diff = Zp - Zbeta_p ;

  if((zp_diff) > 0) {  /* ------------- p-rich */
     delta_m = exp(q1+q2*mass_ratio)* zp_diff ;
  }
  else                /* ------------- n-rich */
  {
    double mass_ratio2 = mass_ratio*mass_ratio;
     delta_m = mass_ratio2*(n1+(n2*mass_ratio2))*zp_diff ;
  }

  Zprob = Zbeta + delta + delta_m + 0.002*A; // Eq (6) from errata
  /* calculate width parameter */
  if ((zp_diff)<0) {           /* n-rich */
     R0 = r[0] * exp(r[3]*zp_diff) ; // Eq (13)
  }
  else  
  {                               /* p-rich */
     R0 = r[0] * exp(r[4]*zp_diff) ; //Eq (14)
  }   
  R = R0 * exp(r[1]*A + r[2]*A*A) ; // Eq (12)

  if (mass_ratio > r[6]) {
    R *= exp(r[5] * sqrt(Ap) * pow(mass_ratio-r[6],3.0)) ; //Eq (15)
  }
 
  if((Zprob-Z) > 0) {
  /*     neutron-rich */
    expo = -R * pow(fabs(Zprob-Z),Un) ;
    fract = exp(expo) * sqrt(R/3.14159) ;
  }
  else {
  /* old V2.1 parameterization */
    double slope = l1 + l2 * pow(A/2.0,0.3) ;
    double z0 = Zprob+( slope*log(10.)/(2.*R)) ;
    expo = -R * pow(fabs(Zprob-Z),Up) ;
    fract   =  exp(expo) * sqrt(R/3.14159) ;
    if( Z > z0 ) {
      expo  = -R * pow(fabs(Zprob-z0),Up) ;
      double  a2    =  exp(expo) * sqrt(R/3.14159) ;
      /* fract =  a2 * exp(slope*(z0-Z)) */
      fract = a2/pow( pow(10,slope),(Z-z0) ) ;
    }
  }

//  fwhm = pow(0.693/R,1.0/Un) + pow(0.693/R,1.0/Up) ;
  norm = 1.0 ;

  /* "brute force" scaling factor for n-rich projectiles only */

  if((zp_diff)<=0) { 
    if ((Zbeta-Z)>(b2+zp_diff)) { 
      norm = pow(10.0, -b1 * fabs(zp_diff) * pow(Zbeta-Z+Zp - Zbeta_p + b2,3)) ;
    }
    else
    {
      norm=1.0 ;
    }
  }

  //double fwhm = pow(0.693/R,1.0/Un) + pow(0.693/R,1.0/Up) ;
  result =  norm * fract * yield_a ; 
  return(result*1000);
}


} // end of epax namespace
