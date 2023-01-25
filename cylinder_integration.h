
#ifndef CYLINDER_INTEGRATION_H
#define CYLINDER_INTEGRATION_H

// #include "integration.h"
#include "cuba_wrapper.h"

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
    
    double operator () (double x[]) // x = {r, phi}
    { 
        return x[0]*(*f_endcap)(x);
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
    
    double operator () (double x[]){ // x = {phi, z}
        return R*(*f_side)(x);
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
                             int mineval, int maxeval, int verbosity,
                             double &result, double &error)
{
    double R = integrand -> R;
    double L = integrand -> L;
    int nregions, neval, fail;
    double zero2[] = {0.0, 0.0};

    (*integrand).surface = top;
    EndcapIntegrand<FuncOnCylinder> top_integrand(R, L, integrand);
    double endcap_uppers[] = {R, 2*M_PI};
    IntegralHC<EndcapIntegrand<FuncOnCylinder>, 2>
        top_integral(&top_integrand, zero2, endcap_uppers);
    double top[1], top_error[1];
    run_cuhre(&top_integral, rtol, atol, verbosity, mineval, maxeval,
              nregions, neval, fail, top, top_error);
    
    (*integrand).surface = bottom;
    EndcapIntegrand<FuncOnCylinder> bottom_integrand(R, L, integrand);
    IntegralHC<EndcapIntegrand<FuncOnCylinder>, 2>
        bottom_integral(&bottom_integrand, zero2, endcap_uppers);
    double bottom[1], bottom_error[1];
    run_cuhre(&bottom_integral, rtol, atol, verbosity, mineval, maxeval,
              nregions, neval, fail, bottom, bottom_error);

    (*integrand).surface = side;
    SideIntegrand<FuncOnCylinder> side_integrand(R, L, integrand);
    double side_uppers[] = {2*M_PI, L};
    IntegralHC<SideIntegrand<FuncOnCylinder>, 2>
        side_integral(&side_integrand, zero2, side_uppers);
    double side[1], side_error[1];
    run_cuhre(&side_integral, rtol, atol, verbosity, mineval, maxeval,
              nregions, neval, fail, side, side_error);
    
    result = top[0] + bottom[0] + side[0];
    error = gsl_hypot3(top_error[0], bottom_error[0], side_error[0]);
}


#endif