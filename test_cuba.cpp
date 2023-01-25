
#include <string>
#include <iostream>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <fmt/format.h>
#include <cuba.h>

#include "cuba_wrapper.h"


double power_integral(double power, double xmin, double xmax){
    double p_plus = power + 1;
    return (pow(xmax, p_plus) - pow(xmin, p_plus))/p_plus;
}


/*
test integrand: x^p1 y^p2, where p1, p2 are parameters
*/
class RaiseToPower2D 
{
public:
    double power0, power1;    
    RaiseToPower2D(double init_power0, double init_power1) : 
        power0(init_power0), power1(init_power1) {}
    
    double operator()(double x[2])
    {
        return pow(x[0], power0)*pow(x[1], power1);
    }
    
    double exact_integral(double x0_min, double x0_max, 
                          double x1_min, double x1_max)
    {
        return (power_integral(power0, x0_min, x0_max)*
                power_integral(power1, x1_min, x1_max)); 
    }
};


/*
test integrand: (x/L) sin(2 pi x y / L H), to be integrated 
over the box (0,0) to (L, H) with L, H > 0. See mathematica.
*/
class MixedSine 
{
public:
    double L, H;    
    MixedSine(double init_L, double init_H) : 
        L(init_L), H(init_H) {}
    
    double operator()(double x[2])
    {
        return (x[0]/L)*gsl_sf_sin(2*M_PI*x[0]*x[1]/(L*H));
    }
    
    /* int_0^H dx1 int_0^L dx2*/
    double exact_integral()
    {
      return L*H/(2*M_PI); 
    }
};


int main(void)
{
    double rtol = 1e-9;
    double atol = 1e-9;
    int mineval = 10;
    int maxeval = 1e9;
    int verbosity = 0;
    int nregions, neval, fail;
    double result[1], error[1];


    std::cout << 
        "\n ----------------- test 2D power law ----------------- " 
             << std::endl;
    RaiseToPower2D cube_square(3, 2);
    double xmin[2] = {0.0, 2.0};
    double xmax[2] = {2.0, 5.0};
    IntegralHC<RaiseToPower2D, 2> 
        integral_cube_square(&cube_square, xmin, xmax);
    run_cuhre(&integral_cube_square, rtol, atol, verbosity, mineval, 
              maxeval, nregions, neval, fail, result, error);
    double actual = cube_square.exact_integral(xmin[0], xmax[0], 
                                               xmin[1], xmax[1]);
    std::cout << fmt::format("numerical = {:f} +- {:e}\n", 
                             result[0], error[0])
              << fmt::format("exact     = {:f} \n", actual)
              << fmt::format("diff      = {:e} +- {:e}\n", 
                             result[0] - actual, error[0])
              << fmt::format("fail? {}\n", fail)
              << fmt::format("nregion = {}\n", nregions) 
              << fmt::format("neval = {}", neval) << std::endl;


    std::cout << 
        "\n ----------------- test mixed sine function ----------------- " 
             << std::endl;
    double L = 3*M_PI*M_PI;
    double R = sqrt(7);
    MixedSine ms(L, R);
    double mins[] = {0.0, 0.0};
    double maxs[] = {L, R};
    IntegralHC<MixedSine, 2> integral_ms(&ms, mins, maxs);
    run_cuhre(&integral_ms, rtol, atol, verbosity, mineval, maxeval,
               nregions, neval, fail, result, error);
    actual = ms.exact_integral();
    std::cout << fmt::format("numerical = {:f} +- {:e}\n", 
                             result[0], error[0])
              << fmt::format("exact     = {:f} \n", actual)
              << fmt::format("diff      = {:e} +- {:e}\n", 
                             result[0] - actual, error[0])
              << fmt::format("fail? {}\n", fail)
              << fmt::format("nregion = {}\n", nregions) 
              << fmt::format("neval = {}", neval) << std::endl;

    return 0;
}
