
#ifndef CYLINDER_INTEGRATION_H
#define CYLINDER_INTEGRATION_H

#include "integration.h"

/*
Map integrand functors (double, double) --> double to new functors of 
the same signature, now including the cylindrical area measure. 
*/
template<typename Func2D>
class EndcapIntegrand {
public:
    double L;
    double R;
    Func2D * f_endcap;
    EndcapIntegrand(double R_init, double L_init, Func2D * f_init) :
      L(L_init), R(R_init), f_endcap(f_init) {}
    double operator () (double r, double phi){
        return r*(*f_endcap)(r, phi);
    }
};

template<typename Func2D>
class SideIntegrand {
public:
    double L;
    double R;
    Func2D * f_side;
    SideIntegrand(double R_init, double L_init, Func2D * f_init) :
      L(L_init), R(R_init), f_side(f_init) {}
    double operator () (double phi, double z){
        return R*(*f_side)(phi, z);
    }
};


/*
Integrate a scalar function over a cylinder. The cylinder has 
coordinates (r, phi, z) with r in [0, R], phi in [0, 2pi), and 
z in [0, L]. The passed FuncOnCyl defines the integrand. This is a
functor (double, double) --> double with parameters 'R' for the
radius, 'L' for the length, and 'surface' which is Surface enum flag 
that selects whether evaluation is to occur on the top, bottom or side 
surface of the cylinder.  The arguments of the functor evaluation 
depend on surface: if stop or bottom, the arguments are assumed to be 
(r, phi) but if surface is side then the arguments are (phi, z).    
*/
enum Surface {top, bottom, side,};

template<typename FuncOnCylinder> 
void integrate_over_cylinder(FuncOnCylinder * integrand, 
                             double atol, double rtol,
                             gsl_integration_method method,
                             double &result, double &error)
{
    double R = integrand -> R;
    double L = integrand -> L;

    (*integrand).surface = top;
    EndcapIntegrand<FuncOnCylinder> top_integrand(R, L, integrand);
    

    
    double top, top_error;
    integrate_over_2d_box(&top_integrand, 
                          0.0, R, 0.0, 2*M_PI, 
                          atol, rtol, method, top, top_error);
    
    (*integrand).surface = bottom;
    EndcapIntegrand<FuncOnCylinder> bottom_integrand(R, L, integrand);
    double bottom, bottom_error;
    integrate_over_2d_box(&bottom_integrand, 
                          0.0, R, 0.0, 2*M_PI, 
                           atol, rtol, method, bottom, bottom_error);

    (*integrand).surface = side;
    SideIntegrand<FuncOnCylinder> side_integrand(R, L, integrand);
    double side, side_error;
    integrate_over_2d_box(&side_integrand,
                          0.0, 2*M_PI, 0.0, L,
                          atol, rtol, method, side, side_error);
    
    result = top + bottom + side;
    error = gsl_hypot3(top_error, bottom_error, side_error);
}


#endif