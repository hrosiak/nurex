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
#ifndef DEFAULT_DENSITIES_H
#define DEFAULT_DENSITIES_H

#include <map>
#include "nurex/Nucleus.h"
#include "nurex/NucID.h"
namespace nurex{


Nucleus get_default_nucleus(int a, int z);

inline Nucleus get_default_nucleus(const char* symbol){
    auto r = nucleus_from_symbol(symbol);
    return get_default_nucleus(r[0],r[1]);
}

}

#endif
