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
#ifndef NUCID_H
#define NUCID_H

#include <map>
#include <string>
#include <algorithm>
#include <array>

namespace nurex{

constexpr inline int nucleus_id(int a, int z){
    return 10000*a + 10*z;
}

constexpr inline std::array<int,2> nucleus_from_id(int id){
    return {id/10000,(id%10000)/10};
}

std::array<int,2> nucleus_from_symbol(std::string symbol);
std::string get_element_symbol(int z);
int get_element_z(std::string el);

inline std::string nucleus_symbol(int a, int z){
    std::string res="";
    if(z>0 && z<119 && a>0){
        res =  std::to_string(a)+get_element_symbol(z);
        }
    return res;
}

}


#endif
