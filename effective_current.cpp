
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_bessel.h>
#include <fmt/format.h>

#include "cylinder_integration.h"
#include "cylinder_modes.h"

#include "effective_current.h"


double propagator(double r, double phi, double z,
                  double r0, double phi0, double z0,
                  double q, PropagatorType prop_type)
{
    double d = sqrt(r*r + r0*r0 - 2.0*r*r0*gsl_sf_cos(phi - phi0) +
                (z - z0)*(z - z0)); // cylindrical coordinates distance
    double arg = q*d; // q must be positive
    switch (prop_type){
    case real:
        return gsl_sf_cos(arg)/d;
    case imaginary:
        return gsl_sf_sin(arg)/d;
    case evanescent:
        if (arg < 100.0){
            return gsl_sf_exp(-arg)/d;
        }
        else { // avoid underflow
            return 0.0;
        }
    }
    abort();
}


PropagatedSurfaceCurrent::PropagatedSurfaceCurrent(
                 double R_init, double L_init,
                 VectorFieldOnCylinder Ki_emitter_init,
                 CylinderFrequency omega_func_init) :
    R(R_init), L(L_init),
    Ki_emitter(Ki_emitter_init),
    omega_func(omega_func_init) {}

/*
Surface current in cylindrical components relative to the
detection point (which is set in the attributes)
*/
double PropagatedSurfaceCurrent::Ki_detector(double x1, double x2, double phi)
{
    switch (component) {
    case r_hat:
        return Ki_emitter(x1, x2, R, L, surface, r_hat)*
                gsl_sf_cos(phi - phi0) +
               Ki_emitter(x1, x2, R, L, surface, phi_hat)*
                gsl_sf_sin(phi - phi0);
    case phi_hat:
        return -Ki_emitter(x1, x2, R, L, surface, r_hat)*
                  gsl_sf_sin(phi - phi0) +
                Ki_emitter(x1, x2, R, L, surface, phi_hat)*
                  gsl_sf_cos(phi - phi0);
    case z_hat:
        return Ki_emitter(x1, x2, R, L, surface, z_hat);
    }
    abort();
}

double PropagatedSurfaceCurrent::operator () (double x1, double x2)
{
    double r, phi, z;
    switch (surface) {
    case top:
        r = x1;
        phi = x2;
        z = L;
        break;
    case bottom:
        r = x1;
        phi = x2;
        z = 0.0;
        break;
    case side:
        r = R;
        phi = x1;
        z = x2;
        break;
    default:
        std::cerr << "invalid unit vector" << std::endl;
        abort();
    }
    double q = wavenumber(omega_func(R, L), m);
    double prop = propagator(r, phi, z, r0, phi0, z0, q, prop_type);
    double Kid = Ki_detector(x1, x2, phi);
    return prop*Kid;
}


EffectiveCurrent::EffectiveCurrent(double R_init, double L_init,
                                   VectorFieldOnCylinder Ki_emitter_init,
                                   CylinderFrequency omega_func_init,
                                   double atol_init, double rtol_init,
                                   gsl_integration_method method_init)
    : mode(R_init, L_init, Ki_emitter_init, omega_func_init),
      atol(atol_init), rtol(rtol_init), method(method_init) {}

void EffectiveCurrent::operator()(double r0, double phi0, double z0,
                                  double m, PropagatorType re_or_im,
                                  CylindricalUnitVector component,
                                  double &result, double &error)
{
    if (r0 < mode.R && 0 < z0 && z0 < mode.L){
        result = NAN;
        error = NAN;
        // std::cout << fmt::format("inside cylinder:\n"
        //                          "  r = {}\n"
        //                          "  z = {}\n\n", r0, z0);
        return;
    }
    mode.r0 = r0;
    mode.phi0 = phi0;
    mode.z0 = z0;
    mode.m = m;
    mode.component = component;
    // select real or imaginary part of j_eff
    double omega = mode.omega_func(mode.R, mode.L);
    if (omega < m){
        if (re_or_im == real){
            mode.prop_type = evanescent;
        } else if (re_or_im == imaginary){
            result = 0.0;
            error = 0.0;
            return;
        }
    } else {
        mode.prop_type = re_or_im;
    }
    integrate_over_cylinder<PropagatedSurfaceCurrent>(&mode, atol, rtol,
                                                      method, result, error);
    result *= -m*m/(4*M_PI);
    error *= -m*m/(4*M_PI);
}
