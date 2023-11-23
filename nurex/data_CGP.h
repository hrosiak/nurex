/*
 *  Copyright(C) 2019, Andrej Prochazka
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

#ifndef DATA_CGP_H
#define DATA_CGP_H
// pairing energies from Cook, Cameron, Gilbert
// taken from Furihata's GEM paper
namespace nurex{
    double pairing_proton[] = { 
    0.00,  5.44,  0.00,  2.76,  0.00,  3.34,  0.00,  2.70,  0.00,  2.50, //1-10
    0.00,  2.46,  0.00,  2.09,  0.00,  1.62,  0.00,  1.62,  0.00,  1.83, //11-20
    0.00,  1.73,  0.00,  1.35,  0.00,  1.54,  0.00,  1.28,  0.26,  0.88, //21-30
    0.19,  1.35, -0.05,  1.52, -0.09,  1.17,  0.04,  1.24,  0.29,  1.09, //31-40
    0.26,  1.17,  0.23,  1.15, -0.08,  1.35,  0.34,  1.05,  0.28,  1.27, //41-50
    0.00,  1.05,  0.00,  1.00,  0.09,  1.20,  0.20,  1.40,  0.93,  1.00, //51-60
   -0.20,  1.19,  0.09,  0.97,  0.00,  0.92,  0.11,  0.68,  0.05,  0.68, //61-70
   -0.22,  0.79,  0.09,  0.69,  0.01,  0.72,  0.00,  0.40,  0.16,  0.73, //71-80
    0.00,  0.46,  0.17,  0.89,  0.00,  0.79,  0.00,  0.89,  0.00,  0.81, //81-90
   -0.06,  0.69, -0.20,  0.71, -0.12,  0.72,  0.00,  0.77 //91-98
    };

    double pairing_neutron[] = { 
    0.00,  5.98,  0.00,  2.77,  0.00,  3.16,  0.00,  3.01,  0.00,  2.50, //1-10
    0.00,  2.67,  0.00,  1.80,  0.00,  1.67,  0.00,  1.86,  0.00,  2.04, //11-20
    0.00,  1.64,  0.00,  1.44,  0.00,  1.54,  0.00,  1.30,  0.00,  1.27, //21-30
    0.00,  1.29,  0.08,  1.41, -0.08,  1.50, -0.05,  2.24, -0.47,  1.43, //31-40
   -0.15,  1.44,  0.06,  1.56,  0.25,  1.57, -0.16,  1.46,  0.00,  0.93, //41-50
    0.01,  0.62, -0.50,  1.42,  0.13,  1.52, -0.65,  0.80, -0.08,  1.29, //51-60
   -0.47,  1.25, -0.44,  0.97,  0.08,  1.65, -0.11,  1.06, -0.46,  1.06, //61-70
    0.22,  1.55, -0.07,  1.37,  0.10,  1.20, -0.27,  0.93, -0.35,  1.19, //71-80
    0.00,  1.05, -0.25,  1.61, -0.21,  0.90, -0.21,  0.74, -0.38,  0.72, //81-90
   -0.34,  0.92, -0.26,  0.94,  0.01,  0.65, -0.36,  0.83,  0.11,  0.67, //91-100
    0.05,  1.00,  0.51,  1.04,  0.33,  0.68, -0.27,  0.81,  0.09,  0.75, //101-110
    0.17,  0.86,  0.14,  1.10, -0.22,  0.84, -0.47,  0.48,  0.02,  0.88, //111- 
    0.24,  0.52,  0.27,  0.41, -0.05,  0.38,  0.15,  0.67,  0.00,  0.61, //121-
    0.00,  0.78,  0.00,  0.67,  0.00,  0.67,  0.00,  0.79,  0.00,  0.60, //131-
    0.04,  0.64, -0.06,  0.45,  0.05,  0.26, -0.22,  0.39,  0.00,  0.39, //141-150
    };

    double cgc_pairing_energy(int A, int Z){
        int N = A - Z;
        double res = 0.0;
        if(Z>0 && Z<=98)res+=pairing_proton[Z-1];
        if(N>0 && N<=150)res+=pairing_neutron[N-1];
        return res;
    }

}
#endif