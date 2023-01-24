
#include <string>
#include <iostream>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <cuba.h>


int test_sine(const int *ndim, const double in[],
              const int *ncomp, double f[], void *userdata)
{
    const double x = in[0];
    const double y = in[1];
    f[0] = gsl_sf_cos(1e4*x)*gsl_sf_exp(y)/sqrt(x*x + 1e-6);
    return 0;
}

int main(void)
{
    int ndim = 2;
    int ncomp = 1;
    integrand_t integrand = &test_sine;
    int empty_params = 0;
    int nvec = 1;
    double epsrel = 1e-9;
    double epsabs = 1e-9;
    int flags = 1;
    int mineval = 10;
    int maxeval = 1e9;
    int key = 13;
    char statefile[] = "";
    int spin_val = -1;
    int nregions, neval, fail;
    double integral[ncomp];
    double error[ncomp];
    double prob[ncomp];

    Cuhre(ndim, ncomp, integrand, &empty_params, nvec,
          epsrel, epsabs, flags, mineval, maxeval, key,
          statefile, &spin_val, &nregions, &neval, &fail,
          integral, error, prob);

    std::cout << std::endl;
    std::cout << integral[0] << std::endl;

    return 0;
}
