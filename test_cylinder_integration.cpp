/*
Test my wrappers for the numerical integration routines of GSL 
*/

#include <string>
#include <iostream>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <fmt/format.h>

#include "cylinder_integration.h"


/**************************************************************************/


/*
test surface area integrand
*/
class Unit2D {
public:
    double R;
    double L;
    Surface surface;

    Unit2D(double R_init, double L_init) : R(R_init), L(L_init) {}
    
    double operator () (double x[]){
        return 1.0;
    }
    
    double exact_surface_area()
    {
        return 2.0*M_PI*R*R + L*2*M_PI*R;
    }
};



/**************************************************************************/


/*
test mixed cosine and exponential integrand
*/
class MixedCosExp {
public:
    double L;
    double R;
    double k;
    Surface surface;

    MixedCosExp(double L_init, double R_init, double k_init) :
      L(L_init), R(R_init), k(k_init) {}

    double bulk(double r, double phi, double z)
    {
        return (z/L)*gsl_sf_sin(z*M_PI/L)*
               gsl_sf_exp(k*r/R)*gsl_sf_cos(phi*z*M_PI/(2*L));
    }

    double exact_suface_integral()
    {
        return 2.0*L*R*gsl_sf_exp(k)*
               gsl_sf_sin(M_PI*M_PI)/((1.0-M_PI*M_PI)*M_PI*M_PI);
    }    

    double operator () (double x[])
    {
        switch (surface){
        case top:
            return bulk(x[0], x[1], L);
        case bottom:
            return bulk(x[0], x[1], 0.0);
        case side:
            return bulk(R, x[0], x[1]);
        }
        abort();
    }
};


/**************************************************************************/


int main (void)
{

    double atol = 1e-12; 
    double rtol = 1e-9; 
    int mineval = 1e3;
    int maxeval = 1e9;
    int verbosity = 0;

    std::cout << 
        "\n ----------------- test surface area ----------------- \n" 
              << std::endl;
    double L = 10.0*M_PI;
    double R = 2.0*pow(2, 0.5);
    Unit2D unit(R, L);
    double result, error;
    integrate_over_cylinder(&unit, atol, rtol, mineval, maxeval, 
                            verbosity, result, error);
    double actual = unit.exact_surface_area();
    std::cout << fmt::format("cuda      = {:f} +- {:e}\n", result, error)
              << fmt::format("exact     = {:f} \n", actual)
              << fmt::format("diff      = {:e} +- {:e}\n", 
                             result - actual, error) 
              << std::endl;


  std::cout << 
    "\n ----------------- test mixed cosine exp ----------------- \n"
            << std::endl;
  double k = -11.0;
  MixedCosExp cos_exp(L, R, k);
  integrate_over_cylinder(&cos_exp, atol, rtol, mineval, maxeval, 
                          verbosity, result, error);
  actual = cos_exp.exact_suface_integral();
  std::cout << fmt::format("\nk = {:f}\n", k)
            << fmt::format("  cuda   = {:e} +- {:e}\n", result, error)
            << fmt::format("  exact     = {:e} \n", actual)
            << fmt::format("  diff      = {:e} +- {:e}\n", 
                           result - actual, error)
            << std::endl;

  k = M_PI*M_PI*2;
  cos_exp.k = k;
  integrate_over_cylinder(&cos_exp, atol, rtol, mineval, maxeval, 
                          verbosity, result, error);
  actual = cos_exp.exact_suface_integral();
  std::cout << fmt::format("\nk = {:f}\n", k)
            << fmt::format("  cuda   = {:e} +- {:e}\n", result, error)
            << fmt::format("  exact     = {:e} \n", actual)
            << fmt::format("  diff      = {:e} +- {:e}\n", 
                           result - actual, error)
            << std::endl;
}