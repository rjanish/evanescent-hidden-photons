
#ifndef CYLINDER_INTEGRATION_H
#define CYLINDER_INTEGRATION_H

#include <cmath>

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
    const double factor = 50.0;

    double R = integrand -> R;
    double L = integrand -> L;

    double omega = (integrand -> omega_func)(integrand -> R, integrand -> L);
    double mass = integrand -> m;
    double range = factor/sqrt(abs(mass*mass - omega*omega));
    double r0 = integrand -> r0;
    double z0 = integrand -> z0;
    double phi0 = integrand -> phi0;
    double r_lower = 0;
    double r_upper = R;
    double z_lower = 0;
    double z_upper = L;
    double phi_lower = 0;
    double phi_upper = 2*M_PI;
    if (mass > omega)
    {
        r_lower = GSL_MAX(0.0, r0 - range);
        r_upper = GSL_MIN(R, r0 + range);
        if (abs(r_lower) > 0)
        {
            double phi_range = range/r_lower;
            if (phi_range < M_PI)
            {
                phi_lower = phi0 - phi_range;
                phi_upper = phi0 + phi_range;
            }
        }
        z_lower = GSL_MAX(0.0, z0 - range);
        z_upper = GSL_MIN(R, z0 + range);
    }

    (*integrand).surface = top;
    EndcapIntegrand<FuncOnCylinder> top_integrand(R, L, integrand);
    double top_result, top_error;
    integrate_over_2d_box(&top_integrand,
                          r_lower, r_upper, phi_lower, phi_upper,
                          atol, rtol, method, top_result, top_error);

    double bottom_result = 0.0;
    double bottom_error = 0.0;
    if (mass < omega || L < range)
    {
        (*integrand).surface = bottom;
        EndcapIntegrand<FuncOnCylinder> bottom_integrand(R, L, integrand);
        integrate_over_2d_box(&bottom_integrand,
                              r_lower, r_upper, phi_lower, phi_upper,
                               atol, rtol, method, bottom_result, bottom_error);
    }

    double side_result = 0.0;
    double side_error = 0.0;
    if (mass < omega || abs(R - r0) < range)
    {
        (*integrand).surface = side;
        SideIntegrand<FuncOnCylinder> side_integrand(R, L, integrand);
        integrate_over_2d_box(&side_integrand,
                              0.0, 2*M_PI, z_lower, z_upper,
                              atol, rtol, method, side_result, side_error);
    }

    result = top_result + bottom_result + side_result;
    error = gsl_hypot3(top_error, bottom_error, side_error);
}


#endif
