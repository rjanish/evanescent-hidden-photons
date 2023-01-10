
#include <iostream>

#include "cylinder_modes.h"
#include "effective_current.h"

#include "overlap.h"


OverlapIntegrand::OverlapIntegrand(double Rs_init, double Ls_init,
                                   VectorFieldOnCylinder Ki_emitter_init,
                                   CylinderFrequency omega_func_init,
                                   double Rd_init, double Ld_init,
                                   double seperation_init,
                                   VectorFieldInCylinder Edetect_init,
                                   PropagatorType re_or_im_init,
                                   double mass_init,
                                   double atol_init, double rtol_init,
                                   gsl_integration_method method_init)
    : j_eff(Rs_init, Ls_init, Ki_emitter_init, omega_func_init,
            atol_init, rtol_init, method_init),
      Rd(Rd_init), Ld(Ld_init), seperation(seperation_init),
      Edetect(Edetect_init), mass(mass_init), re_or_im(re_or_im_init),
      atol(atol_init), rtol(rtol_init), method(method_init) {}


double OverlapIntegrand::operator()(double r, double z)
{
    double j_dot_E = 0.0;
    const double phi = 0.0; // assume azimuthal symmetry
    const CylindricalUnitVector directions[3] = {r_hat, phi_hat, z_hat};
    for(auto &spatial_component : directions)
    {
        double Ji, error;
        j_eff(r, phi, j_eff.mode.L + seperation + z,
              mass, re_or_im, spatial_component, Ji, error);
        j_dot_E += Ji*Edetect(r, phi, z, Rd, Ld, spatial_component);
    }
    return 2*M_PI*r*j_dot_E;
}


Overlap::Overlap(double Rs_init, double Ls_init,
                 VectorFieldOnCylinder Ki_emitter_init,
                 CylinderFrequency omega_func_init,
                 double Rd_init, double Ld_init, double seperation_init,
                 VectorFieldInCylinder Edetect_init,
                 double atol_init, double rtol_init,
                 gsl_integration_method method_init)
    : integrand(Rs_init, Ls_init, Ki_emitter_init, omega_func_init,
                Rd_init, Ld_init, seperation_init, Edetect_init,
                real, 0.0, // placeholders
                atol_init, rtol_init, method_init),
      atol(atol_init), rtol(rtol_init), method(method_init) {}


double Overlap::operator()(double mass)
{
    integrand.mass = mass;
    const PropagatorType re_or_im[2] = {real, imaginary};
    double overlap_sq = 0.0;
    for(auto &complex_part : re_or_im){
        integrand.re_or_im = complex_part;
        double result, error;
        integrate_over_2d_box(&integrand,
                              0.0, integrand.Rd,
                              0.0, integrand.Ld,
                              atol, rtol, method, result, error);
        overlap_sq += result*result;
    }
    return sqrt(overlap_sq);
}

