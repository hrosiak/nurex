/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 */
//#include <limits>
#include "spline.h"
namespace nurex {

namespace vector{
tridiagonal_matrix::tridiagonal_matrix(unsigned int dim)
{
    a.resize(dim);
    d.resize(dim);
    c.resize(dim);
}

// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & tridiagonal_matrix::operator () (unsigned int i, unsigned int j)
{
    int k=j-i;
    if(k == -1)return c[i];
    else if(k==0) return d[i];
    else return a[i];
}

double tridiagonal_matrix::operator () (unsigned int i, unsigned int j) const
{
    int k=j-i;
    if(k==-1)return c[i];
    else if(k==0) return d[i];
    else if(k==1)return a[i];
    else return 0.0;
}


std::vector<double> tridiagonal_matrix::trig_solve(const std::vector<double>& b) const
{
    std::vector<double> x(this->dim());
    if(b.size() != d.size() || d[0] == 0.0){return x;}
    std::vector<double> g(this->dim());
    x[0] = b[0]/d[0];
    double bet = d[0];
    for(std::size_t j=1, max=d.size();j<max;j++){
        g[j] = c[j-1]/bet;
        bet = d[j] - (a[j]*g[j]);
        if(bet == 0.0){
	    x.clear();
            return x;
        }
        x[j] = (b[j]-a[j]*x[j-1])/bet;
    }
    for(int j=dim()-2;j>=0;j--){
        x[j] -= g[j+1]*x[j+1];
    }
    return x;
}

  }

cspline_vector::cspline_vector(const std::vector<double>& x, const std::vector<double>& y, bool boundary_second_deriv){
    assert(x.size()==y.size());
    assert(x.size()>2);
    m_x=x;
    m_y=y;
    auto n=x.size();

    vector::tridiagonal_matrix A(n);
    std::vector<double>  rhs(n);
    for(std::size_t i=1; i<n-1; i++) {
        A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
        A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
        A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
        assert((x[i+1]-x[i])!=0);
        assert((x[i]-x[i-1])!=0);
        rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    }
        // boundary conditions
    if(boundary_second_deriv) {
        // 2*b[0] = f''
        A(0,0)=2.0;
        A(0,1)=0.0;
        rhs[0]=0.0; // 0.0 is value of derivative
        A(n-1,n-1)=2.0;
        A(n-1,n-2)=0.0;
        rhs[n-1]=0.0; // 0.0 is value of derivative
       } else {
        // c[0] = f', needs to be re-expressed in terms of b:
        // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
        A(0,0)=2.0*(x[1]-x[0]);
        A(0,1)=1.0*(x[1]-x[0]);
        rhs[0]=3.0*((y[1]-y[0])/(x[1]-x[0])-0.0); // 0.0 is deriv value

        // c[n-1] = f', needs to be re-expressed in terms of b:
        // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
        // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
        A(n-1,n-1)=2.0*(x[n-1]-x[n-2]);
        A(n-1,n-2)=1.0*(x[n-1]-x[n-2]);
        rhs[n-1]=3.0*(0.0-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }


        m_b=A.trig_solve(rhs);

        // calculate parameters a[] and c[] based on b[]
        m_a.resize(n);
        m_c.resize(n);
        for(int i=0; i<n-1; i++) {
            m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
            m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
                   - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
        }


    // for left extrapolation coefficients
    //s.m_b0 = (m_force_linear_extrapolation==false) ? s.m_b[0] : 0.0;
    m_b0 =  0.0;
    m_c0 = m_c[0];

    double h=x[n-1]-x[n-2];
    m_a[n-1]=0.0;
    m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
    m_b[n-1]=0.0;
};


double cspline_vector::evaluate(double x) const
{
    size_t n=m_x.size();
    // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
    std::vector<double>::const_iterator it;
    it=std::lower_bound(m_x.begin(),m_x.end(),x);
    int idx=std::max( int(it-m_x.begin())-1, 0);

    double h=x-m_x[idx];
    double interpol;
    if(x<m_x[0]) {
        // extrapolation to the left
        interpol=(m_b0*h + m_c0)*h + m_y[0];
    } else if(x>m_x[n-1]) {
        // extrapolation to the right
        interpol=(m_b[n-1]*h + m_c[n-1])*h + m_y[n-1];
    } else {
        // interpolation
        interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
    }

    return interpol;
}

}


