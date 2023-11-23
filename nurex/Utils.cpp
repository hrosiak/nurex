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
#include "nurex/Utils.h"
namespace nurex{

double finite_range(double r, double b)
{
    double res;
    double b2 = b * b;
    res = 0.5 * exp(-0.5 * r * r / b2) / (nurex::PI * b2);
    return res;
}

double finite_range2(double r, double b)
{
    double res;
    double b2 = b * b;
    res = 0.5 * exp(-0.5 * r / b2) / (nurex::PI * b2);
    return res;
}


std::vector<std::pair<int,double>> get_discrete_vector_by_area(std::pair<double,double> par){
    using std::abs;
    std::vector<std::pair<int,double>> res;
    double sum = 0.0;
    double mean = par.first;
    double sigma = par.second;
    double norm = std::sqrt(2)*sigma;
    if(sigma<0.1){
        res.push_back({mean,1.0});
        return res;
    }
    for(int x=static_cast<int>(std::floor(mean-2*sigma)); x<=static_cast<int>(std::ceil(mean+2*sigma)); x = x+1 ){
        double p = x-mean;

     /*
        const double d = 0.5;
        double arg = p;
        double argp = (arg>=0.0)?arg+d:arg-d;
        if(x==0)arg=(arg>=0.0)?arg-d:arg+d;
 */
        double arg = p - 0.5;
        double argp = p + 0.5;

        double a;
        if(arg*argp>0){
            a = 0.5*abs( std::erf(abs(argp)/norm) - std::erf(abs(arg)/norm) );
        }
        else {
            a = 0.5*(std::erf(abs(argp)/norm) + std::erf(abs(arg)/norm));
        }
//        a = gaussian(p,mean,sigma);
        //printf("%d %lf %lf -> %lf\n",x,arg,argp,a);
                sum+=a;
        res.push_back({x,a});
    }
    if(sum<0.99 || sum>1.01){
        for(auto &e:res){
            e.second = e.second/sum;
        }
    }
    return res;
}



} //end of sigmar namespace
