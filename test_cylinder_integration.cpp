/*
Test my wrappers for the numerical integration routines of GSL 
*/

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>

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
    double operator () (double x, double y){
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

    double operator () (double x1, double x2)
    {
        switch (surface){
        case top:
            return bulk(x1, x2, L);
        case bottom:
            return bulk(x1, x2, 0.0);
        case side:
            return bulk(R, x1, x2);
        }
        abort();
    }
};


/**************************************************************************/


int main (void)
{

  printf("\n ----------------- test surface area ----------------- \n");
  double L = 10.0*M_PI;
  double R = 2.0*pow(2, 0.5);
  Unit2D unit(R, L);

  double actual = unit.exact_surface_area();
  double result, error;
  double atol = 1e-12; 
  double rtol = 1e-9; 
  integrate_over_cylinder<Unit2D>(&unit, atol, rtol, romberg, result, error);
  printf("romberg   = %f +- %e\n", result, error);
  printf("exact     = %f \n", actual);
  printf("diff      = %e +- %e\n\n", result - actual, error);
  integrate_over_cylinder<Unit2D>(&unit, atol, rtol, adaptive_singular, 
                                  result, error);
  printf("qags      = %f +- %e\n", result, error);
  printf("exact     = %f \n", actual);
  printf("diff      = %e +- %e\n", result - actual, error);


  printf("\n\n ----------------- test mixed cosine exp ----------------- \n");
  double k = -11.0;
  MixedCosExp cos_exp(L, R, k);
  integrate_over_cylinder<MixedCosExp>(&cos_exp, atol, rtol, 
                                       romberg, result, error);
  actual = cos_exp.exact_suface_integral();
  printf("\nk = %f\n", k);
  printf("\tromberg   = %e +- %e\n", result, error);
  printf("\texact     = %e \n", actual);
  printf("\tdiff      = %e +- %e\n\n", result - actual, error);
  integrate_over_cylinder<MixedCosExp>(&cos_exp, atol, rtol, 
                                       adaptive_singular, result, error);
  printf("\tqags      = %e +- %e\n", result, error);
  printf("\texact     = %e \n", actual);
  printf("\tdiff      = %e +- %e\n", result - actual, error);

  k = M_PI*M_PI*2;
  cos_exp.k = k;
  integrate_over_cylinder<MixedCosExp>(&cos_exp, atol, rtol, romberg, result, error);
  actual = cos_exp.exact_suface_integral();
  printf("\nk = %f\n", k);
  printf("\tromberg   = %e +- %e\n", result, error);
  printf("\texact     = %e \n", actual);
  printf("\tdiff      = %e +- %e\n\n", result - actual, error);
  integrate_over_cylinder<MixedCosExp>(&cos_exp, atol, rtol, 
                                       adaptive_singular, result, error);
  printf("\tqags      = %e +- %e\n", result, error);
  printf("\texact     = %e \n", actual);
  printf("\tdiff      = %e +- %e\n", result - actual, error);
}