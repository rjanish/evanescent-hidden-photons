#ifndef CUBA_WRAPPER_H
#define CUBA_WRAPPER_H


#include <cuba.h>


/*
Define an integral over a ndim hypercube, with the integrand a vector 
of ncomp components. integrand is a callable object with signiture 
  void integrand(double input[ndim], double output[ncomp]). 
xmin and xmax are (ndim,) arrays giving the boundaries of the 
integration region.  
*/
template<typename Func, int ndim_in>
class IntegralHC
{
public:
    int ndim = ndim_in; 
    Func * integrand; 
    double xmin[ndim_in];
    double xmax[ndim_in];
    double xlength[ndim_in];
    double scale = 1.0;

    IntegralHC(Func * integrand_init, double xmin_init[], double xmax_init[]) 
    : integrand(integrand_init)
    {
        for (int index=0; index < ndim_in; ++index) 
        {
            xmin[index] = xmin_init[index];
            xmax[index] = xmax_init[index];
            xlength[index] = xmax[index] - xmin[index];
            scale *= xlength[index];
        }
    }

    void integrand_from_unitHC(const double in_unitHC[], double out[])
    {
        double in_physical[ndim_in]; 
        for (int index=0; index < ndim_in; ++index)
        {
            in_physical[index] = 
                xmin[index] + in_unitHC[index]*xlength[index];
        }
        out[0] = (*integrand)(in_physical);
    }
};


/*
This is a function on the unit hypercube that follows cuda's integrand_t
signiture.  It is passed a pointer to an IntegeralHC object which is
used to evaluate the integrand and scale the argument from the unit 
hypercube into the space expected by IntegralHC.integrand. 
*/
template<typename Func, int ndim>
int cuba_dressed_integrand(const int *ndim_dummy,  // unused here
                           const double input_unitHC[],
                           const int *ncomp_dummy, // unused here
                           double out[], void * ptr_to_IntegralHC)
{
    IntegralHC<Func, ndim> * integral_definition = 
        static_cast<IntegralHC<Func, ndim> *>(ptr_to_IntegralHC);
    (integral_definition -> integrand_from_unitHC)(input_unitHC, out);
    return 0;
};


template<typename Func, int ndim>
void run_cuhre(IntegralHC<Func, ndim> * integral_definition, 
               double rtol, double atol, int verbosity, int mineval, 
               int maxeval, int &nregions, int &neval, int &fail,
               double result[], double error[])
{
    int ncomp = 1;  // only run one component at a time for now
    int nvec = 1;           // evaluates integrand one at a time
    int key = 7;            // use default order for quadrature rule
    char statefile[] = "";  // do not save the integration state to file
    int spin_val = -1;      // handle threading automatically
    double prob[ncomp];
    integrand_t integrand_wrapper = &cuba_dressed_integrand<Func, ndim>;
    double scale = integral_definition -> scale;
    Cuhre(ndim, ncomp, integrand_wrapper, integral_definition, 
          nvec, rtol, atol/scale, verbosity, mineval, maxeval, key,
          statefile, &spin_val, &nregions, &neval, &fail,
          result, error, prob);
    for (int index=0; index < ncomp; ++index)
    {
        result[index] *= scale;
        error[index]  *= scale;
    }
};


#endif
