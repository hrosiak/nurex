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
#include "nurex/NNCrossSection.h"
#include <fstream>

namespace nurex{

///////////// ------ NNCrossSectionFit ----- /////////////////

double NNCrossSectionFit::sigma_np(double x) const{
    double result = 0.0;
    // This will be fits to the np cross sections
    if(x<0.00882){
      result = 20360.;
    }

    else if(x>=0.00882 && x<0.0505){
      result = -160030.+168000.*pow(x,-0.015033);
    }

    else if(x>=0.0505 && x<0.2){
      result = -163380.+166200.*pow(x,-0.02494);
    }

    else if(x>=0.2 && x<0.38){
      result = (-163380.+166200.*pow(x,-0.02494))*(0.38-x)/0.18+(-503.08+4761.3*pow(x,-0.49562))*(x-0.2)/0.18;
      }

    //else if(x>=0.38 && x<2.){
    else if(x>=0.38 && x<1.5479){
      result = -503.08+4761.3*pow(x,-0.49562);
    }

    else if(x>=1.5479 && x<2.){
      result = (-503.08+4761.3*pow(x,-0.49562))*(2.-x)/0.4521+(POWER5(8.3738,-0.63495,0.14901,-0.1317,0.033962,-0.0034031,x))*(x-1.5479)/0.4521;
      }

    else if(x>=2 && x<34.){
      result = POWER5(8.3738,-0.63495,0.14901,-0.1317,0.033962,-0.0034031,x);
    }

    else if(x>=34. && x<598.){
      result = POWER5(7.8594,1.291,-0.95288,0.13433,-0.0057926,0.000092646,x);
    }

    else if(x>=598 && x<700.){
      result = (POWER5(7.8594,1.291,-0.95288,0.13433,-0.0057926,0.000092646,x))*(700.-x)/102.+(POWER5(-112.96,25.957,1.0306,-0.086902,-0.09889,0.0090234,x))*(x-598.)/102.;
      }

    //    else if(x>=600. && x<1000.){
    else if(x>=700. && x<981.18){
      result = POWER5(-112.96,25.957,1.0306,-0.086902,-0.09889,0.0090234,x);
    }

    else if(x>=981.18 && x<2500.){
      result = POWER5(-7.3768,0.43052,0.4123,0.037464,-0.018571,0.0011638,x);
    }

    else{
      result = POWER5(-7.3768,0.43052,0.4123,0.037464,-0.018571,0.0011638,2500.);
    }

    return result;
}
double NNCrossSectionFit::sigma_pp(double x) const{
    double result=0;

    // This will be fits to the pp cross sections
    if(x<1.){
      result = 0.;
    }
    else if(x>=1. && x<1.5){
      result = (-503.08+4761.3*(pow(x,-0.49562)))/2.917;
        }

    else if(x>=1.5 && x<2.5){
      result = ((-503.08+4761.3*(pow(x,-0.49562)))/2.917)*(2.5-x)/1.+(POWER5(8.3738,-0.63495,0.14901,-0.1317,0.033962,-0.0034031,x)/2.917)*(x-1.5)/1.;
      }

    else if(x>=2.5 && x<11.05){
      result = POWER5(8.3738,-0.63495,0.14901,-0.1317,0.033962,-0.0034031,x)/2.917;
        }

    else if(x>=11.05 && x<50.){
      result = POWER5(13.257,-5.832,1.2461,0.077779,-0.074157,0.0078546,x);
        }

    else if(x>=50. && x<70.){
      result = (POWER5(13.257,-5.832,1.2461,0.077779,-0.074157,0.0078546,x))*(70.-x)/20.+(POWER5(-43.793,31.3,-2.2222,-2.308,0.5649,-0.037926,x))*(x-50.)/20.;
      }

    else if(x>=70. && x<150.37){
      result = POWER5(-43.793,31.3,-2.2222,-2.308,0.5649,-0.037926,x);
        }

    else if(x>=150.37 && x<260.87){
      result = -0.0076434*x+25.5503;
      }

    else if(x>=260.87 && x<411.970){
      result = POWER5(-8.5816,3.3171,0.29625,-0.056291,-0.024141,0.0031779,x);
        }
    else if(x>=411.970 && x<600.){
      result = EXPO(20.826,0.49428,-0.0057621,x);
        }
    else if(x>=600. && x<680.){
      result = (EXPO(20.826,0.49428,-0.0057621,x))*(680.-x)/80.+(EXPO(47.614,-1.2511e5,0.015477,x))*(x-600.)/80.;
      }
    else if(x>=680. && x<1000.){
      result = EXPO(47.614,-1.2511e5,0.015477,x);
        }
    else if(x>=1000. && x<3000.){
      result = POLY(45.847,0.0052698,-4.1863e-6,6.8537e-10,x);
        }
    else{
      result = POLY(45.847,0.0052698,-4.1863e-6,6.8537e-10,3000.);
        }
    return result;
}

double NNCrossSectionFit::np(double E){
    #ifdef USE_THREADS
    std::lock_guard<std::mutex> lock(mutex);
    #endif
    if(std::fabs(E-energy_np)<0.0001)return cs_np;

    cs_np = sigma_np(E);
    energy_np = E;

    return cs_np;
  }

  double NNCrossSectionFit::pp(double E){
    #ifdef USE_THREADS
    std::lock_guard<std::mutex> lock(mutex);
    #endif
    if(std::fabs(E-energy_pp)<0.0001)return cs_pp;
    cs_pp = sigma_pp(E);
    energy_pp = E;
    return cs_pp;
  }

  double NNCrossSectionFit::EXPO(double AA, double BB, double CC, double energy)noexcept{
    double result = AA + BB*std::exp(-CC*energy);
    return result;
  }

  double NNCrossSectionFit::POLY(double AA, double BB, double CC, double DD, double energy)noexcept{
    double result = AA + BB*energy + CC*energy*energy + DD*energy*energy*energy;
    return result;
  }

  double NNCrossSectionFit::POWER5(double AA, double BB, double CC, double DD, double EE, double FF, double energy)noexcept{
    double loge = log(energy);
    double loge2 = loge*loge;
    double result = std::exp(AA + BB*loge
          + CC*loge2
          + DD*loge2*loge
          + EE*loge2*loge2
          + FF*loge2*loge2*loge);
    return result;

  }

#ifndef NO_FILESYSTEM
///////////// ------ NNCrossSectionFile ----- /////////////////
NNCrossSectionFile::NNCrossSectionFile(const char *np_filename, const char *pp_filename){
    std::string line;
    std::ifstream fpn,fpp;
    fpn.open(np_filename);
    fpp.open(pp_filename);

    if(!fpn.is_open() || !fpp.is_open()){
        throw std::invalid_argument("Failed to open the file");
    }

    npfile = np_filename;
    ppfile = pp_filename;

    double _x,_y, _e;
    std::vector<double> x,y;

    while(fpn>>_x && fpn>>_y && fpn>>_e){
        if(x.size()>0 && _x<=x.back()){
            printf("not increasing energy:  %lf, prev.:%lf \n",_x,x.back());
        }
        if(_y<0){
            printf("cross section less than 0, energy:  %lf, cc: %lf \n",_x,_y);
        }
        x.push_back(_x);
        y.push_back(_y);
    }
    fpn.close();
    //printf("loaded %u data points\n",x.size());
    ipn = new nurex::Interpolator(x,y);
    x.clear();
    y.clear();

    while(fpp>>_x && fpp>>_y && fpp>>_e){
        if(x.size()>0 && _x<=x.back()){
            printf("not increasing energy:  %lf, prev.:%lf \n",_x,x.back());
        }
        if(_y<0){
            printf("cross section less than 0, energy:  %lf, cc: %lf \n",_x,_y);
        }
        x.push_back(_x);
        y.push_back(_y);
    }
    //printf("loaded %u data points\n",x.size());
    ipp = new nurex::Interpolator(x,y);
    fpp.close();
    x.clear();
    y.clear();
}

NNCrossSectionFile::~NNCrossSectionFile(){
    delete ipn;
    delete ipp;
    ipn=nullptr;
    ipp=nullptr;
}

double NNCrossSectionFile::np(double energy){
    return (*ipn)(energy);
}

double NNCrossSectionFile::pp(double energy){
    return (*ipp)(energy);
}
#endif

template class FermiMotion<NNCrossSectionFit>;
template class FermiMotionD<NNCrossSectionFit>;

} //end of nurex namespace
