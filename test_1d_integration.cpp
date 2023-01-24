/*
Test my wrappers for the numerical integration routines of GSL
*/

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_trig.h>

#include "integration.h"


/*
test integrand: x^p, p is a parameter
*/
class RaiseToPower {
  public:
    double power;
    RaiseToPower(double ap) : power(ap) {}
    double operator () (double base){
      return pow(base, power);
    }
};

/*
test integrand: sin(x) or cos(x) based on selection parameter
*/
class Trig {
  public:
    char s_or_c;
    Trig(char selector) : s_or_c(selector) {}
    double operator () (double x){
      double output;
      if (s_or_c == 's'){
        output = gsl_sf_sin(x);
      }
      else if (s_or_c == 'c'){
        output = gsl_sf_cos(x);
      }
      else{
        abort();
      }
      return output;
    }
};




int main (void)
{

  // test power function

  RaiseToPower cube(3);
  gsl_function integrand_quad = make_gsl_function_from_functor(&cube);

  printf("\n ----- power func ----- \n");
  printf("4^3 ?= %f\n", cube(4.0));

  double result, error;
  size_t neval;

  gsl_integration_qng(&integrand_quad, 0, 1, 1e-12, 1e-7, &result, &error, &neval);

  printf("int_0^1 dx x^3 ?= %f\n", result);

  // test trig function
  printf("\n ----- trig funcs ----- \n");

  Trig trig('c');
  gsl_function integrand_trig = make_gsl_function_from_functor(&trig);

  printf("cos(0) ?= %f\n", trig(0.0));
  printf("cos(pi/2) ?= %f\n", trig(M_PI/2.0));

  gsl_integration_qng(&integrand_trig, 0, M_PI, 1e-12, 1e-7, &result, &error, &neval);
  printf("int_0^pi dx cos(x) ?= %f\n", result);
  gsl_integration_qng(&integrand_trig, 0, M_PI/2.0, 1e-12, 1e-7, &result, &error, &neval);
  printf("int_0^pi/2 dx cos(x) ?= %f\n", result);

  trig.s_or_c = 's';  // switch to sin
  printf("\n");
  printf("sin(0) ?= %f\n", trig(0.0));
  printf("sin(pi/2) ?= %f\n", trig(M_PI/2.0));

  gsl_integration_qng(&integrand_trig, 0, M_PI, 1e-12, 1e-7, &result, &error, &neval);
  printf("int_0^pi dx sin(x) ?= %f\n", result);
  gsl_integration_qng(&integrand_trig, 0, M_PI/2.0, 1e-12, 1e-7, &result, &error, &neval);
  printf("int_0^pi/2 dx sin(x) ?= %f\n", result);

  return 0;
}
