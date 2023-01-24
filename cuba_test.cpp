
#include <string>
#include <iostream>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>

#include "cuba_wrapper.h"


class TestSine
{
public:
    double k;
    double d;
    int ndim = 2;
    int ncomp = 1;
    TestSine(double k_init, double d_init) : k(k_init), d(d_init) {};
    int operator()(double in[2], double out[1])
    {
        auto x = in[0];
        auto y = in[1];
        out[0] = gsl_sf_cos(k*x)*gsl_sf_exp(y)/sqrt(x*x + d*d);
        return 0;
    }
};



int main(void)
{
    // // int ndim = 2;
    // // int ncomp = 1;
    // // integrand_t integrand = &test_sine;
    // // int empty_params = 0;
    // // int nvec = 1;
    // double rtol = 1e-9;
    // double atol = 1e-9;
    // // int flags = 0b100;
    // int mineval = 10;
    // int maxeval = 1e9;
    // int verbosity = 1;
    // // int key = 13;
    // // char statefile[] = "";
    // // int spin_val = -1;
    // int nregions, neval, fail;
    double integral[1];
    // double error[1];

    TestSine fast_sine(1e4, 1e-3);
    std::cout << integral[0] << std::endl;
    double xtest[]={0,0};
    fast_sine(xtest, integral);
    std::cout << integral[0] << std::endl;

    // Cuhre(ndim, ncomp, integrand, &empty_params, nvec,
    //       epsrel, epsabs, flags, mineval, maxeval, key,
    //       statefile, &spin_val, &nregions, &neval, &fail,
    //       integral, error, prob);

    // run_cuhre<TestSine>(&fast_sine, rtol, atol,
    //                     verbosity, mineval, maxeval,
    //                     &nregions, &neval, &fail,
    //                     integral, error);

    // std::cout << std::endl;
    // std::cout << integral[0] << std::endl;
    // std::cout << error[0] << std::endl;

    return 0;
}
