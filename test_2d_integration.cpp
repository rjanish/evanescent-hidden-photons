/*
Test my wrappers for the numerical integration routines of GSL 
*/

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_trig.h>

#include "integration.h"


double power_integral(double power, double xmin, double xmax){
    double p_plus = power + 1;
    return (pow(xmax, p_plus) - pow(xmin, p_plus))/p_plus;
}

/*
test integrand: x^p1 y^p2, where p1, p2 are parameters
*/
class RaiseToPower2D {
  public:
    double power0, power1;
    RaiseToPower2D(double init_power0, double init_power1) : 
      power0(init_power0), power1(init_power1) {}
    double operator () (double x0, double x1){
      return pow(x0, power0)*pow(x1, power1);
    }
    double exact_integral(double x0_min, double x0_max, 
                          double x1_min, double x1_max){
      return (power_integral(power0, x0_min, x0_max)*
               power_integral(power1, x1_min, x1_max)); 
    }
};


/*
test integrand: (x/L) sin(2 pi x y / L H), to be integrated 
over the box (0,0) to (L, H) with L, H > 0. See mathematica.
*/
class MixedSin {
  public:
    double L, H;
    MixedSin(double init_L, double init_H) : 
      L(init_L), H(init_H) {}
    double operator () (double x, double y){
      return (x/L)*gsl_sf_sin(2*M_PI*x*y/(L*H));
    }
    double exact_integral_from_origin()
    {
      return L*H/(2*M_PI); 
    }
};


int main (void)
{

  printf("\n ----------------- test power functor ----------------- \n");
  RaiseToPower2D cube_square(3, 2);
  double x, y;
  x = 2;
  y = 3;
  printf("(%.1f)^3 (%.1f)^2 = %f\n", x, y, cube_square(x, y));
  x = 10;
  y = 4;
  printf("(%.1f)^3 (%.1f)^2 = %f\n", x, y, cube_square(x, y));

  printf("\n --- test argument fixer --- \n");
  FixFirstArgument<RaiseToPower2D> square(&cube_square);
  square.x0 = 1.0;
  printf("(%.1f)^2 = %f\n", y, square(y));

  printf("\n --- test integration: romberg --- \n");
  double x0_min = M_PI;
  double x0_max = 11*M_PI;
  double x1_min = -2.0;
  double x1_max = pow(13, 0.5);
  double atol = 1e-12; 
  double rtol = 1e-9; 
  double result, error;
  integrate_over_2d_box<RaiseToPower2D>(&cube_square, 
                                        x0_min, x0_max,
                                        x1_min, x1_max,
                                        atol, rtol, romberg,
                                        result, error);
  double actual = cube_square.exact_integral(x0_min, x0_max, x1_min, x1_max);
  printf("numerical = %f +- %e\n", result, error);
  printf("exact     = %f \n", actual);
  printf("diff      = %e +- %e\n", result - actual, error);

  printf("\n --- test integration: adaptive_singular --- \n");
  integrate_over_2d_box<RaiseToPower2D>(&cube_square, 
                                        x0_min, x0_max,
                                        x1_min, x1_max,
                                        atol, rtol, adaptive_singular,
                                        result, error);
  printf("numerical = %f +- %e\n", result, error);
  printf("exact     = %f \n", actual);
  printf("diff      = %e +- %e\n", result - actual, error);


  printf("\n");
  printf("\n ----------------- test mixed sine functor ----------------- \n");
  double L = 15;
  double H = 5;
  MixedSin ms(L, H);
  x = 2;
  y = 3;
  printf("ms(%.1f, %.1f) = %f\n", x, y, ms(x, y));
  x = 4.9;
  y = 13.3;
  printf("ms(%.1f, %.1f) = %f\n", x, y, ms(x, y));

  printf("\n --- test argument fixer --- \n");
  FixFirstArgument<MixedSin> ms_y(&ms);
  ms_y.x0 = 2;
  y = 3;
  printf("ms(%.1f, %.1f) = %f\n", ms_y.x0, y, ms_y(y));

  printf("\n --- test integration: romberg --- \n");
  x0_min = 0.0;
  x0_max = L;
  x1_min = 0.0;
  x1_max = H;
  atol = 1e-12; 
  rtol = 1e-9; 
  integrate_over_2d_box<MixedSin>(&ms, 
                                  x0_min, x0_max,
                                  x1_min, x1_max,
                                  atol, rtol, romberg,
                                  result, error);
  actual = ms.exact_integral_from_origin();
  printf("numerical = %f +- %e\n", result, error);
  printf("exact     = %f \n", actual);
  printf("diff      = %e +- %e\n", result - actual, error);

  printf("\n --- test integration: adaptive_singular --- \n");
  integrate_over_2d_box<MixedSin>(&ms, 
                                  x0_min, x0_max,
                                  x1_min, x1_max,
                                  atol, rtol, adaptive_singular,
                                  result, error);
  printf("numerical = %f +- %e\n", result, error);
  printf("exact     = %f \n", actual);
  printf("diff      = %e +- %e\n", result - actual, error);
}