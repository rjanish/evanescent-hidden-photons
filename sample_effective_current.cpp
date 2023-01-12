
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

    // fixed parameters
    const gsl_integration_method method = adaptive_singular;
    gsl_set_error_handler_off();
    const size_t num_input_params = 15;
    const CylindricalUnitVector directions[3] = {r_hat, phi_hat, z_hat};
    const PropagatorType re_or_im[2] = {real, imaginary};

    // process each param file
    for (int file_number = 1; file_number < argc; ++file_number){
        // read parameter file
          // z here is in 'detector coordinates', i.e. z measured
          // from the top of the source cylinder: z_detector = L + z_source
          // phi in the param file should be given in units of 2pi, but
          // in the output file it will be given in radians
        std::string param_filename = argv[file_number];
        std::ifstream param_file(param_filename);
        std::map<std::string, std::string> input_params;
        read_param_entries(param_file, input_params, num_input_params);
        std::string mode_name = input_params["mode"];
        auto r_N     = std::stoul(input_params["r_N"]);
        auto r_min   = std::stod(input_params["r_min"]);
        auto r_max   = std::stod(input_params["r_max"]);
        auto phi_N   = std::stoul(input_params["phi_N"]);
        auto phi_min = std::stod(input_params["phi_min"])*2*M_PI;
        auto phi_max = std::stod(input_params["phi_max"])*2*M_PI;
        auto z_N     = std::stoul(input_params["z_N"]);
        auto z_min   = std::stod(input_params["z_min"]);
        auto z_max   = std::stod(input_params["z_max"]);
        auto radius  = std::stod(input_params["radius"]);
        auto length  = std::stod(input_params["length"]);
        auto atol    = std::stod(input_params["atol"]);
        auto rtol    = std::stod(input_params["rtol"]);
        auto mass    = std::stod(input_params["mass"]);
        std::cout << fmt::format("{}...", param_filename) << std::endl;

        // construct sampling grid
        double r_grid[r_N];
        linspace(r_min, r_max, r_N, r_grid);
        double phi_grid[phi_N];
        linspace(phi_min, phi_max, phi_N, phi_grid);
        double z_grid[z_N];
        linspace(z_min, z_max, z_N, z_grid);

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
        std::ofstream output_file(fmt::format("output-{}", param_filename));
        write_map(output_file, input_params);
        output_file << '\n';
        output_file << fmt::format("r  phi  z  "
                                   "Re[j_r]  error_Re[j_r]  "
                                   "Im[j_r]  error_Im[j_r]  "
                                   "Re[j_phi]  error_Re[j_phi]  "
                                   "Im[j_phi]  error_Im[j_phi]  "
                                   "Re[j_z]  error_Re[j_z]  "
                                   "Im[j_z]  error_Im[j_z]  \n{}\n\n",
                                   std::string(78, '-'));

        // sample effective current
        double result, error;
        for (auto &r : r_grid){
            for (auto &phi : phi_grid){
                for (auto &z : z_grid){
                    output_file << fmt::format("{}  {}  {}  ", r, phi, z);
                    for(auto &spatial_component : directions){
                        for(auto &complex_part : re_or_im){
                            J(r, phi, length + z, mass,
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
    return 0;
}
