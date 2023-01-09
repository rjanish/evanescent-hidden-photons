
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>

#include "cylinder_modes.h"
#include "effective_current.h"


int main()
{
    double R = 1.0;
    double L = 3.0;
    double atol = 1e-12; 
    double rtol = 1e-9; 
    gsl_integration_method method = romberg;
    EffectiveCurrent J(R, L, &Ki_cylinder_TE011, &angular_frequency_TE011, 
                       atol, rtol, method);
    double r0 = 0.5;
    double phi0 = 0.0;
    double z0 = L + 1e-4;
    CylindricalUnitVector component = phi_hat;

    double m_min = 0.05;
    double m_max = 50;
    int N_m = 8; // number of points between endpoints, ie: total - 2 
    int N_m_total = N_m + 2;
    double d_logm = (gsl_sf_log(m_max) - gsl_sf_log(m_min))/(N_m_total - 1);

    double result, error;
    for (int index=0; index < N_m_total; ++index) {
        double m = m_min*gsl_sf_exp(d_logm*index);

        J(r0, phi0, z0, m, real, component, result, error);
        printf("Re[j(%f)] = %e +- %e\n", m, result, error);

        J(r0, phi0, z0, m, imaginary, component, result, error);
        printf("Im[j(%f)] = %e +- %e\n\n", m, result, error);
    }

    return 0;
}