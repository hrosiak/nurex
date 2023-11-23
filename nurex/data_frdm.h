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

#ifndef FRDM_H
#define FRDM_H

namespace frdm{
constexpr inline int nucleus_id(int a, int z){
    return 10000*a + 10*z;
}

struct data{
    int nid;
    double b2;
    double b3;
    double b4;
    double EFLmic;
};

double get_b2(int nid);
double get_b3(int nid);
double get_b4(int nid);
double get_EFLmic(int nid);
const data& get_data(int nid);

inline const data& get_data(int A, int Z){
    return get_data(nucleus_id(A,Z));
}
    

inline double get_b2(int A, int Z){
    return get_data(nucleus_id(A,Z)).b2;
}

inline double get_b3(int A, int Z){
    return get_data(nucleus_id(A,Z)).b3;
}

inline double get_b4(int A, int Z){
    return get_data(nucleus_id(A,Z)).b4;
}

inline double get_EFLmic(int A, int Z){
    return get_data(nucleus_id(A,Z)).EFLmic;
}

} //end of namespace
#endif
