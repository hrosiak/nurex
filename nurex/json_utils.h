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
#ifndef JSON_UTILS_H
#define JSON_UTILS_H
#include <vector>
#include "nurex/Density.h"
#include "nurex/GlauberModel.h"
namespace nurex{

GlauberModelType json_model(const std::string& string);
std::vector<double> json_energies(const std::string& string);
DensityType json_density(const std::string& string);
EvaporationParameters json_evaporation_parameters(const std::string& string);
}
#endif
