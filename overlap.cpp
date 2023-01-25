
#include <iostream>

#include <gsl/gsl_sf_exp.h>

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
                                   int mineval_init, int maxeval_init,
                                   int verbosity_init)
    : mode(Rs_init, Ls_init, Ki_emitter_init, omega_func_init),
      Rd(Rd_init), Ld(Ld_init), seperation(seperation_init),
      Edetect(Edetect_init), mass(mass_init), re_or_im(re_or_im_init),
      atol(atol_init), rtol(rtol_init), mineval(mineval_init), 
      maxeval(maxeval_init), verbosity(verbosity_init) {}


double OverlapIntegrand::operator()(double x[])
{
    mode.prop_type = re_or_im;
    mode.m = mass;
    
    double x_mode[] = {x[0], x[1]};
    double r = x[2];
    mode.r0 = r;
    double phi = 0.0; // assume azimuthal symmetry
    mode.phi0 = phi;
    double z = x[3];
    mode.z0 = z;
    
    const CylindricalUnitVector directions[3] = {r_hat, phi_hat, z_hat};
   
    double j_dot_E = 0.0;
    for(auto &spatial_component : directions)
    {
        mode.component = spatial_component;
        j_dot_E += mode(x_mode)*Edetect(r, phi, z, Rd, Ld, spatial_component);
    }
    return r*j_dot_E;
}


Overlap::Overlap(double Rs_init, double Ls_init,
                 VectorFieldOnCylinder Ki_emitter_init,
                 CylinderFrequency omega_func_init,
                 double Rd_init, double Ld_init, double seperation_init,
                 VectorFieldInCylinder Edetect_init,
                 double atol_init, double rtol_init,
                 int mineval_init, int maxeval_init, int verbosity_init)
    : integrand(Rs_init, Ls_init, Ki_emitter_init, omega_func_init,
                Rd_init, Ld_init, seperation_init, Edetect_init,
                real, 0.0, // placeholders
                atol_init, rtol_init, 
                mineval_init, maxeval_init, verbosity_init),
      atol(atol_init), rtol(rtol_init), mineval(mineval_init), 
      maxeval(maxeval_init), verbosity(verbosity_init) {}


double Overlap::operator()(double mass)
{
    integrand.mass = mass;
    const PropagatorType re_or_im[2] = {real, imaginary};
    double overlap_sq = 0.0;
    for(auto &complex_part : re_or_im){
        integrand.re_or_im = complex_part;
        double result, error;
        integrate_over_cylinderXplane(&integrand, atol, rtol, mineval, 
                                      maxeval, verbosity, method, 
                                      result, error);
        overlap_sq += result*result; 
            // integral over r and z, assuming integrand is independent of phi
    }
    auto volume = M_PI*integrand.Rd*integrand.Rd*integrand.Ld;
    auto omega = integrand.j_eff.mode.omega_func(integrand.Rd, integrand.Ld);
    auto scale = mass*gsl_sf_exp(mass*integrand.seperation)/sqrt(omega*volume);
     // This normalization make the overlap factor unitless and inpedendent of m 
     // for large m, and removed the overall exponential gap scaling. See paper.
    return -scale*0.5*M_PI*m*m*sqrt(overlap_sq); // 0.5 = 2pi/4pi 

