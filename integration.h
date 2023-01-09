
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>


/*
Convert any functor f: double --> double to a gsl_function function 
object. The gsl_function will be given no extra parameters, these 
must be set via the functor's attributes. 
*/
template<typename Func>
double functor_to_function(double x, void * rawptr)
{
  Func * f = static_cast<Func *>(rawptr);
  return (*f)(x);
};

template<typename Func>
gsl_function make_gsl_function_from_functor(Func * f){
  gsl_function new_gsl_function;
  new_gsl_function.function = &functor_to_function<Func>;
  new_gsl_function.params = f;
  return new_gsl_function;
};


/*
Convert any functor f: (double, double) --> double to a single
argument functor f_fixed double --> double whose one argument is 
the old functor's second argument. This just holds a pointer to the 
original functor and the omitted first arguments as a parameter. 
*/
template<typename Func2D>
class FixFirstArgument {
private:
    Func2D * f;
public:
    double x0 = NAN;
    FixFirstArgument(Func2D * f_init) : f(f_init) {}
    double operator () (double x1){
        return (*f)(x0, x1);
    }
};


/*
Function to quickly switch gsl integration methods
*/
enum gsl_integration_method {adaptive_singular, romberg};

void gsl_integrate(gsl_function * , double, double, double, double, 
                   gsl_integration_method, double &, double &);

/*
This is functor which computes the innermost integral of a 2D box integral,
which is assumed to be over the second argument of the integrand:
  II(x0) = int_x1_min^x1_max dx1 integrand(x0, x1) 
The full integral is given by integrating this over the first argument:
  int_x0_min^x0_max dx0 int_x1_min^x1_max dx1 integrand(x0, x1) = 
  int_x0_min^x0_max dx0 II(x0) 
The integrand is a functor integrand: (double, double) --> double.
*/
template<typename Func2D>
class InnermostIntegration {
private:
    double x1_min;
    double x1_max;
    double atol;
    double rtol;
    gsl_integration_method method;
    FixFirstArgument<Func2D> integrand_with_fixed_first_arg;
    gsl_function gsl_integrand;
public: 
    InnermostIntegration(Func2D * integrand_init,
                         double x1_min_init, double x1_max_init,
                         double atol_init, double rtol_init,
                         gsl_integration_method method_init) :
      x1_min(x1_min_init), x1_max(x1_max_init), 
      atol(atol_init), rtol(rtol_init), method(method_init),
      integrand_with_fixed_first_arg(
        FixFirstArgument<Func2D>(integrand_init)),
      gsl_integrand(
        make_gsl_function_from_functor(&integrand_with_fixed_first_arg)) {}    
    double operator () (double x0) {
        integrand_with_fixed_first_arg.x0 = x0;
        double result, error;
        gsl_integrate(&gsl_integrand , x1_min, x1_max, atol, rtol, 
                      method, result, error);
        return result;
    }
};


/*
Compute the integral of a functor integrand: (double, double) --> double
over a 2D box [x0_min, x0_max]x[x1_min, x1_max] where x0 is the first 
argument to integrand and x1 the second.
*/
template<typename Func2D> 
void integrate_over_2d_box(Func2D * integrand, 
                           double x0_min, double x0_max,
                           double x1_min, double x1_max,
                           double atol, double rtol,
                           gsl_integration_method method,
                           double &result, double &error)
{   
    InnermostIntegration<Func2D> inner_integral(integrand, x1_min, x1_max, 
                                                atol, rtol, method);
    gsl_function gsl_integrand = 
        make_gsl_function_from_functor(&inner_integral);
    gsl_integrate(&gsl_integrand , x0_min, x0_max, atol, rtol, 
                  method, result, error);
};


#endif