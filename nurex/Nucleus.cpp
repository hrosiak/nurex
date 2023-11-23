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
#include "nurex/Nucleus.h"
#include <stdexcept>
#include "nurex/DefaultDensities.h"

namespace nurex{


Nucleus::Nucleus(int _A, int _Z, DensityType const &_density_p, DensityType const &_density_n):a(_A),z(_Z){
    density_p = _density_p;
    density_n = _density_n;
    if(density_p.Norm()!=z){
        density_p.Normalize(z);
    }
    if(density_n.Norm()!=N()){
        density_n.Normalize(N());
    }
}

/**
 * Nucleus copy constructor
 */
Nucleus::Nucleus(const Nucleus &nuc){
    if(!nuc){
        throw std::invalid_argument("Nucleus class not propertly initialized");
    }
    density_p = nuc.density_p;
    density_n = nuc.density_n;
    a = nuc.A();
    z = nuc.Z();
}

/**
 *  Nucleus move constructor
 */
Nucleus::Nucleus(Nucleus &&nuc){
    if(!nuc){
        throw std::invalid_argument("Nucleus class not propertly initialized");
    }
    a = nuc.A();
    z = nuc.Z();
    density_p = std::move(nuc.density_p);
    density_n = std::move(nuc.density_n);
}

Nucleus::Nucleus(const char* sym){
    Nucleus temp = get_default_nucleus(sym);
    if(!temp){
        throw std::invalid_argument("Nucleus class not propertly initialized");
    }
    a = temp.A();
    z = temp.Z();
    density_p = std::move(temp.density_p);
    density_n = std::move(temp.density_n);
}

Nucleus& Nucleus::operator=(Nucleus&& src){
    density_p = std::move(src.density_p);
    density_n = std::move(src.density_n);
    a = src.A();
    z = src.Z();
    return *this;
}


/**
 * bool operator
 */
Nucleus::operator bool() const{
    if(a==-1 || !density_p || !density_n)
        return false;
    else
        return true;
}

////////////// Compound //////////////////
Compound::Compound(std::initializer_list<atom>list){
    std::initializer_list<atom>::iterator it;
    atoms.reserve(list.size());
    for ( it=list.begin(); it!=list.end(); ++it){
        add_element(*it);
    }
}

void Compound::add_element(const atom &a){
    atoms.push_back(a);
    molar_mass += a.stn*a.nucleus.A();
}

}
