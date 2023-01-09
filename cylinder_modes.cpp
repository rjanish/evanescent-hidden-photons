
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>

#include "cylinder_integration.h"

#include "cylinder_modes.h"


/********************* TM010 ****************************************/


double angular_frequency_TM010(double R, double L)
{
    const double x01 = gsl_sf_bessel_zero_J0(1);
    return x01/R;
}


double Ki_cylinder_TM010(double x1, double x2, double R, double L,
                         Surface surface, CylindricalUnitVector component)
{
    const double x01 = gsl_sf_bessel_zero_J0(1);
    switch (surface) {
    case top:
        // x1 is r, x2 is phi
        switch (component) {
        case r_hat:
            return gsl_sf_bessel_J1(x01*x1/R);
        case phi_hat:
            return 0.0;
        case z_hat:
            return 0.0;
        }
        abort();
    case bottom:
        // x1 is r, x2 is phi
        switch (component) {
        case r_hat:
            return -gsl_sf_bessel_J1(x01*x1/R);
        case phi_hat:
            return 0.0;
        case z_hat:
            return 0.0;
        }
        abort();
    case side:
        // x1 is phi, x2 is z
        switch (component) {
        case r_hat:
            return 0.0;
        case phi_hat:
            return 0.0;
        case z_hat:
            return -gsl_sf_bessel_J1(x01);
        }
        abort();
    }
    abort();
}


double Ei_cylinder_TM010(double r, double phi, double z,
                         double R, double L, CylindricalUnitVector component)
    {
        switch (component) {
        case r_hat:
            return 0.0;
        case phi_hat:
            return 0.0;
        case z_hat:
            const double x01 = gsl_sf_bessel_zero_J0(1);
            return -gsl_sf_bessel_J0(x01*r/R);
        }
        abort()
    }


/********************* TE011 ****************************************/


double angular_frequency_TE011(double R, double L)
{
    const double xprime01 = gsl_sf_bessel_zero_J1(1);  // x'_01 = x_11
    return gsl_hypot(xprime01/R, M_PI/L);
}


double Ki_cylinder_TE011(double x1, double x2, double R, double L,
                         Surface surface, CylindricalUnitVector component)
{
    const double xprime01 = gsl_sf_bessel_zero_J1(1);  // x'_01 = x_11
    switch (surface) {
    case top:
        // x1 is r, x2 is phi
        switch (component) {
        case r_hat:
            return 0.0;
        case phi_hat:
            return M_PI*R*gsl_sf_bessel_J1(xprime01*x1/R)/(xprime01*L);
        case z_hat:
            return 0.0;
        }
        abort();
    case bottom:
        // x1 is r, x2 is phi
        switch (component) {
        case r_hat:
            return 0.0;
        case phi_hat:
            return M_PI*R*gsl_sf_bessel_J1(xprime01*x1/R)/(xprime01*L);
        case z_hat:
            return 0.0;
        }
        abort();
    case side:
        // x1 is phi, x2 is z
        switch (component) {
        case r_hat:
            return 0.0;
        case phi_hat:
            return -gsl_sf_bessel_J0(xprime01)*gsl_sf_sin(M_PI*x2/L);
        case z_hat:
            return 0.0;
        }
        abort();
    }
    abort();
}


double Ei_cylinder_TE011(double r, double phi, double z,
                         double R, double L, CylindricalUnitVector component)
    {
        switch (component) {
        case r_hat:
            return 0.0;
        case phi_hat:
            const double omega = angular_frequency_TE011(R, L);
            const double xprime01 = gsl_sf_bessel_zero_J1(1);
            return (omega*R/xprime01) *
                   gsl_sf_bessel_J1(xprime01*r/R) *
                   gsl_sf_sin(M_PI*z/L);
                   // I have dropped a factor of -i here, I think that is fine
        case z_hat:
            return 0.0;
        }
        abort()
    }
