
#include <string>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <fmt/format.h>

#include "cylinder_modes.h"
#include "effective_current.h"
#include "datafiles.h"


int main()
{
    auto R = 1.0;
    auto L = 3.0;
    auto atol = 1e-12; 
    auto rtol = 1e-9; 
    gsl_integration_method method = adaptive_singular;
    EffectiveCurrent J(R, L, &Ki_cylinder_TE011, &angular_frequency_TE011, 
                       atol, rtol, method);
    auto r0 = 0.5;
    auto phi0 = 0.0;
    auto z0 = L + 1e-4;
    CylindricalUnitVector component = phi_hat;

    auto m_min = 0.05;
    auto m_max = 50;
    size_t N_m = 10;

    double result, error;
    for (size_t index=0; index < N_m; ++index) {
        auto m = log_step(m_min, m_max, N_m, index);

        J(r0, phi0, z0, m, real, component, result, error);
        std::cout << fmt::format("Re[j({:f})] = {:e} +- {:e}\n", m, result, error);

        J(r0, phi0, z0, m, imaginary, component, result, error);
        std::cout << fmt::format("Im[j({:f})] = {:e} +- {:e}\n\n", m, result, error);
    }

    return 0;
}