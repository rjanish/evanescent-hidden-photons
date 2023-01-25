#ifndef CUBA_WRAPPER_H
#define CUBA_WRAPPER_H


#include <cuba.h>


/*
Convert any functor f(in[ndim], out[ncomp]) to a function
with a cuda integrand signature. All extra parameters should be
set in the functor's attributes, including ndim and ncomp.
*/
template<typename Func>
int cuba_dressed_integrand(const int *ndim, const double x[],
                   const int *ncomp, double f[],
                   void * ptr_to_integrand)
{
    Func * integrand = static_cast<Func *>(ptr_to_integrand);
    (*integrand)(x, f);
    return;
};


template<typename Func>
void run_cuhre(Func * integrand, double rtol, double atol,
               int verbosity, int mineval, int maxeval,
               int &nregions, int &neval, int &fail,
               double integral[], double error[])
{
    int ndim = integrand -> ndim;
    int ncomp = integrand -> ncomp;
    int nvec = 1;          // evaluates integrand one at a time
    int flags = verbosity;  // 0 to 3
    int key = 13;          // order of quadrature rule
    char statefile[] = ""; // do not save the integration state to file
    int spin_val = -1;     // handle threading automatically
    double prob[ncomp];

    Cuhre(ndim, ncomp, &cuba_dressed_integrand, integrand, nvec,
          rtol, atol, flags, mineval, maxeval, key,
          statefile, &spin_val, &nregions, &neval, &fail,
          integral, error, prob);
};


#endif
