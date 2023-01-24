
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include <gsl/gsl_errno.h>
#include <fmt/format.h>

#include "effective_current.h"
#include "datafiles.h"


void linspace(double start, double end, size_t N, double samples[])
{
    if (N == 1){
        samples[0] = start;
        return;
    } else if (N > 1){
        double step = (end - start)/(N - 1);
        for (size_t index = 0; index < N; ++index){
            samples[index] = start + step*index;
        }
        return;
    } else {
        std::cerr << "invalid sample number" << std::endl;
        abort();
    }
}


int main(int argc, char* argv[]){
    gsl_set_error_handler_off();
    const gsl_integration_method method = adaptive_singular;
    const CylindricalUnitVector directions[3] = {r_hat, phi_hat, z_hat};
    const PropagatorType re_or_im[2] = {real, imaginary};

    // process each param file
    for (int file_number = 1; file_number < argc; ++file_number){
        std::string param_filename = argv[file_number];
        std::cout << fmt::format("{}...", param_filename) << std::endl;
        // read parameter file
          // z here is in 'detector coordinates', i.e. z measured
          // from the top of the source cylinder: z_detector = L + z_source
          // phi in the param file should be given in units of 2pi, but
          // in the output file it will be given in radians
        std::ifstream param_file(param_filename);
        std::map<std::string, std::string> input_params;
        std::string comment("#");
        read_param_entries(param_file, input_params, comment);
        std::string mode_name = input_params["mode"];
        auto m_N     = std::stoul(input_params["m_N"]);
        auto m_start   = std::stod(input_params["m_start"]);
        auto m_end   = std::stod(input_params["m_end"]);
        auto r_N     = std::stoul(input_params["r_N"]);
        auto r_start   = std::stod(input_params["r_start"]);
        auto r_end   = std::stod(input_params["r_end"]);
        auto phi_N   = std::stoul(input_params["phi_N"]);
        auto phi_start = std::stod(input_params["phi_start"])*2*M_PI;
        auto phi_end = std::stod(input_params["phi_end"])*2*M_PI;
        auto z_N     = std::stoul(input_params["z_N"]);
        auto z_start   = std::stod(input_params["z_start"]);
        auto z_end   = std::stod(input_params["z_end"]);
        auto radius  = std::stod(input_params["radius"]);
        auto length  = std::stod(input_params["length"]);
        auto atol    = std::stod(input_params["atol"]);
        auto rtol    = std::stod(input_params["rtol"]);

        // construct effective current functor
        VectorFieldOnCylinder Ki_emitter;
        CylinderFrequency omega_func;
        if (mode_name == "TE011"){
            Ki_emitter = &Ki_cylinder_TE011;
            omega_func = &angular_frequency_TE011;
        } else if (mode_name == "TM010"){
            Ki_emitter = &Ki_cylinder_TM010;
            omega_func = &angular_frequency_TM010;
        } else {
            std::cerr << fmt::format("invalid mode: {}\n", mode_name);
            abort();
        }
        EffectiveCurrent J(radius, length, Ki_emitter, omega_func,
                           atol, rtol, method);

        // prep output file
        std::ofstream output_file(fmt::format("{}.out", param_filename));
        input_params["omega"] = fmt::format("{:0.6e}", omega_func(radius, length));
        write_map(output_file, input_params);
        output_file << '\n';
        output_file << fmt::format("m    r    phi    z    "
                                   "Re[j_r]    error_Re[j_r]    "
                                   "Im[j_r]    error_Im[j_r]    "
                                   "Re[j_phi]    error_Re[j_phi]    "
                                   "Im[j_phi]    error_Im[j_phi]    "
                                   "Re[j_z]    error_Re[j_z]    "
                                   "Im[j_z]    error_Im[j_z]  \n{}\n\n",
                                   std::string(78, '-'));

        // sample effective current
        double result, error;
        for (size_t i_m=0; i_m < m_N; ++i_m){
            double m = log_step(m_start, m_end, m_N, i_m);
            for (size_t i_r=0; i_r < r_N; ++i_r){
                double r = linear_step(r_start, r_end, r_N, i_r);
                for (size_t i_phi=0; i_phi < phi_N; ++i_phi){
                    double phi = linear_step(phi_start, phi_end, phi_N, i_phi);
                    for (size_t i_z=0; i_z < z_N; ++i_z){
                        double z = linear_step(z_start, z_end, z_N, i_z);
                        output_file << fmt::format("{}  {}  {}  {}  ",
                                                   m, r, phi, z);
                        for(auto &spatial_component : directions){
                            for(auto &complex_part : re_or_im){
                                J(r, phi, length + z, m,
                                  complex_part, spatial_component,
                                  result, error);
                                output_file << fmt::format("{}  {}  ",
                                                           result, error);
                            }
                        }
                        output_file << '\n';
                    }
                }
            }
        }
    }
    return 0;
}
