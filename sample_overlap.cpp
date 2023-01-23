
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include <fmt/format.h>
#include <gsl/gsl_errno.h>

#include "overlap.h"
#include "effective_current.h"
#include "cylinder_modes.h"
#include "datafiles.h"


int main(int argc, char* argv[]){

    // fixed parameters
    const gsl_integration_method method = adaptive_singular;
    gsl_set_error_handler_off();

    // process each param file
    for (int file_number = 1; file_number < argc; ++file_number){
        std::string param_filename = argv[file_number];
        std::cout << fmt::format("{}...", param_filename) << std::endl;
        // read parameter file
        std::ifstream param_file(param_filename);
        std::map<std::string, std::string> input_params;
        std::string comment("#");
        read_param_entries(param_file, input_params, comment);
        std::string source_mode = input_params["mode_s"];
        std::string detect_mode = input_params["mode_d"];
        auto Rs      = std::stod(input_params["Rs"]);
        auto Rd      = std::stod(input_params["Rd"]);
        auto Ls      = std::stod(input_params["Ls"]);
        auto Ld      = std::stod(input_params["Ld"]);
        auto atol    = std::stod(input_params["atol"]);
        auto rtol    = std::stod(input_params["rtol"]);
        auto sep     = std::stod(input_params["sep"]);
        auto m_start = std::stod(input_params["m_start"]);
        auto m_end   = std::stod(input_params["m_end"]);
        auto m_N     = std::stoul(input_params["m_N"]);

        // construct overlap function
        VectorFieldOnCylinder Ki_emitter;
        CylinderFrequency omega_func;
        if (source_mode == "TE011"){
            Ki_emitter = &Ki_cylinder_TE011;
            omega_func = &angular_frequency_TE011;
        } else if (source_mode == "TM010"){
            Ki_emitter = &Ki_cylinder_TM010;
            omega_func = &angular_frequency_TM010;
        } else {
            std::cerr << fmt::format("invalid mode: {}\n", source_mode);
            abort();
        }

        VectorFieldInCylinder Edetect;
        if (source_mode == "TE011"){
            Edetect = &Ei_cylinder_TE011;
        } else if (source_mode == "TM010"){
            Edetect = &Ei_cylinder_TM010;
        } else {
            std::cerr << fmt::format("invalid mode: {}\n", detect_mode);
            abort();
        }

        Overlap overlap(Rs, Ls, Ki_emitter, omega_func,
                        Rd, Ld, sep, Edetect, atol, rtol, method);

        // prep output file
        std::ofstream output_file(fmt::format("{}.out", param_filename));
        input_params["omega"] = omega_func(Rs, Ls);
        write_map(output_file, input_params);
        output_file << '\n';
        output_file << fmt::format("m    overlap \n{}\n\n",
                                   std::string(78, '-'));

        double mass, result;
        for (unsigned long index=0; index < m_N; ++index){
            mass = log_step(m_start, m_end, m_N, index);
            std::cout << fmt::format("{}/{}:  m = {}", index+1, m_N, mass)
                      << std::endl;
            result = overlap(mass);
            output_file << fmt::format("{}  {}\n", mass, result);
        }
    }
    return 0;
}
