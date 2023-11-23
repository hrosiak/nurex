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
#ifndef GHQ_INTEGRATOR_H
#define GHQ_INTEGRATOR_H
#include <array>

namespace integrators{

constexpr double sqrt2 = 1.41421356237309504880;
constexpr double sqrt2_inv = 1.0/sqrt2;
constexpr double pi = 3.14159265358979323846;
constexpr double sqrt_pi_inv = 0.56418958354;

template<int N>
struct GH_data{
};

template<int order>
class GaussHermiteIntegration{
public:
    template<typename F>
    double integrate(F& f, double mean=0.0, double sigma=sqrt2_inv, bool norm=false) const;
    template<typename F>
    double operator()(F& f, double mean=0.0, double sigma = sqrt2_inv, bool norm=false) const {return integrate(f, mean, sigma, norm);}
    double w(int i) const {
        return GH_data<order/2>::w()[i];
        }
    double x(int i) const {
        return GH_data<order/2>::x()[i];
    }
    int n() const {return order;}
    std::array<double,order> get_points(double mean=0.0, double sigma=sqrt2_inv)const;
};

template<int order>
template<typename F>
double GaussHermiteIntegration<order>::integrate(F& f, double mean, double sigma, bool norm) const{
    double res=0.0;
    double p = sqrt2*sigma;
    for(int i=0;i<order/2;i++){
        res += w(i) * (f(p*x(i) + mean) + f(-p*x(i) + mean));
    }
    //return res*sqrt_pi_inv;
    return (norm)?res*sqrt_pi_inv:res*sqrt2*sigma;
}

template<int order>
std::array<double,order> GaussHermiteIntegration<order>::get_points(double mean, double sigma)const{
    std::array<double,order> points;
    double p = sqrt2*sigma;
    int num = (order/2);
    for(int i=0;i< num;i++){
        points[num-i-1] = -p*x(i) + mean;
        points[num+i] = p*x(i) + mean;
    }
    return points;
}


template<int order>
class GaussHermiteIntegration2D{
public:
    template<typename F>
    double integrate(F& f, double mean1=0.0, double sigma1=sqrt2_inv, double mean2=0.0, double sigma2=sqrt2_inv) const;
    template<typename F>
    double operator()(F& f, double mean1=0.0, double sigma1=sqrt2_inv, double mean2=0.0, double sigma2=sqrt2_inv) const {return integrate(f, mean1, sigma1, mean2, sigma2);}
private:
    GaussHermiteIntegration<order> integrator;
};

template<int order>
template<typename F>
double GaussHermiteIntegration2D<order>::integrate(F& f, double mean1, double sigma1, double mean2, double sigma2) const{
    double res=0.0;
    double p = sqrt2*sigma1;
    double r = sqrt2*sigma2;
    double xx, yy, sum;
    for(int i=0;i<order/2;i++){
        xx = (p*integrator.x(i));
        for(int j=0;j<order/2;j++){
            yy = (r*integrator.x(j));
            sum = (f(xx + mean1, yy + mean2) + f(xx + mean1, -yy+mean2));
            sum += (f(-xx + mean1, yy + mean2) + f(-xx + mean1, -yy+mean2));
            res += integrator.w(i) * integrator.w(j) * sum;
        }
    }
    return 2.0*res*sigma1*sigma2;
    return res;
}

/// weights and abscissas
// order = 4
template<>
struct GH_data<2>{
    static std::array<double,2> const & x(){
        static const std::array<double,2> _x = {0.5246476232752903178841,1.650680123885784555883};
        return _x;
    }
    static std::array<double,2> const & w(){
        static const std::array<double,2> _w = {0.8049140900055128365061, 0.08131283544724517714304};
        return _w;
        }
};


// order = 6
template<>
struct GH_data<3>{
    static std::array<double,3> const & x(){
        static const std::array<double,3> _x = {0.436077411927616508679, 1.335849074013696949715,2.350604973674492222834};
        return _x;
    }
    static std::array<double,3> const & w(){
        static const std::array<double,3> _w = {0.72462959522439252409, 0.1570673203228566439163, 0.004530009905508845640858};
        return _w;
    }
};

// order = 8
template<>
struct GH_data<4>{
    static std::array<double,4> const & x(){
        static const std::array<double,4> _x = {0.3811869902073221168547, 1.157193712446780194721, 1.981656756695842925855, 2.930637420257244019224};
        return _x;
    }
    static std::array<double,4> const & w(){
        static const std::array<double,4> _w = {0.6611470125582412910304, 0.2078023258148918795433, 0.0170779830074134754562, 1.996040722113676192061E-4};
        return _w;
    }
};


// order = 10
template<>
struct GH_data<5>{
    static std::array<double,5> const & x(){
        static const std::array<double,5> _x = {0.3429013272237046087892, 1.036610829789513654178,1.756683649299881773451,2.532731674232789796409,3.436159118837737603327};
        return _x;
    }
    static std::array<double,5> const & w(){
        static const std::array<double,5> _w = {0.6108626337353257987836, 0.2401386110823146864165,0.03387439445548106313617,0.001343645746781232692202,7.64043285523262062916E-6};
        return _w;
    }
};

// order = 16
template<>
struct GH_data<8>{
    static std::array<double,8> const & x(){
        static const std::array<double,8> _x = {0.27348104613815, 0.82295144914466, 1.3802585391989, 1.9517879909163, 2.5462021578475, 3.17699916198, 3.8694479048601, 4.6887389393058};
        return _x;
    }
    static std::array<double,8> const & w(){
        static const std::array<double,8> _w = {0.50792947901661, 0.28064745852853, 0.083810041398986, 0.01288031153551, 9.322840086242E-4, 2.711860092538E-5,2.320980844865E-7, 2.654807474011E-10 };
        return _w;
    }
};

//oder = 20
template<>
struct GH_data<10>{
    static std::array<double,10> const & x(){
        static const std::array<double,10> _x = {
            0.2453407083009012499038,
            0.7374737285453943587056,
            1.234076215395323007886,
            1.738537712116586206781,
            2.254974002089275523082,
            2.78880605842813048053,
            3.347854567383216326915,
            3.944764040115625210376,
            4.603682449550744273078,
            5.387480890011232862017};
        return _x;
    }
    static std::array<double,10> const & w(){
        static const std::array<double,10> _w = {
            0.46224366960061008965,
            0.28667550536283412972,
            0.1090172060200233200138,
            0.0248105208874636108822,
            0.00324377334223786183218,
            2.283386360163539672572E-4,
            7.8025564785320636941E-6,
            1.086069370769281694E-7,
            4.39934099227318055363E-10,
            2.22939364553415129252E-13,
            };
        return _w;
    }
};

//oder = 32
template<>
struct GH_data<16>{
    static std::array<double,16> const & x(){
        static const std::array<double,16> _x = {
            0.1948407415693993267087,
            0.584978765435932448467,
            0.976500463589682838485,
            1.370376410952871838162,
            1.767654109463201604628,
            2.169499183606112173306,
            2.577249537732317454031,
            2.992490825002374206285,
            3.417167492818570735874,
            3.853755485471444643888,
            4.305547953351198445263,
            4.777164503502596393036,
            5.275550986515880127819,
            5.812225949515913832766,
            6.409498149269660412174,
            7.1258139098307275728
            };
        return _x;
    }
    
    static std::array<double,16> const & w(){
        static const std::array<double,16> _w = {
            0.375238352592802392867,
            0.277458142302529898138,
            0.1512697340766424825752,
            0.06045813095591261418659,
            0.01755342883157343030344,
            0.00365489032665442807913,
            5.3626836552797204597E-4,
            5.416584061819982558E-5,
            3.65058512956237605737E-6,
            1.57416779254559402927E-7,
            4.09883216477089661824E-9,
            5.9332914633966386145E-11,
            4.21501021132644757297E-13,
            1.197344017092848665829E-15,
            9.2317365365182922335E-19,
            7.3106764273841623933E-23
            };
        return _w;
    }
};


}//end of namespace
#endif
