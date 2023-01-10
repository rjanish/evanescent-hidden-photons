
#include <string>
#include <iostream>

#include <fmt/format.h>
#include <gsl/gsl_errno.h>

#include "overlap.h"
#include "effective_current.h"
#include "cylinder_modes.h"


int main(int argc, char* argv[]){


    const gsl_integration_method method = adaptive_singular;
    gsl_set_error_handler_off();
    double atol = 1e-9;
    double rtol = 1e-6;

    double Rs = 1.0;
    double Ls = 3.0;
    VectorFieldOnCylinder Ki_emitter = &Ki_cylinder_TE011;
    CylinderFrequency omega_func = &angular_frequency_TE011;
    double Rd = 1.0;
    double Ld = 3.0;
    VectorFieldInCylinder Edetect = &Ei_cylinder_TE011;
    double seperation = 1e-4;

    Overlap overlap(Rs, Ls, Ki_emitter, omega_func,
                    Rd, Ld, seperation, Edetect, atol, rtol, method);


    double mass = 0.01;
    double result = overlap(mass);

    std::cout << fmt::format("overlap({}) = {}", mass, result);

    return 0;
}
