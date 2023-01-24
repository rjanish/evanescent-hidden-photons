
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "integration.h"


void gsl_integrate(gsl_function * gsl_func , double x_min, double x_max,
                   double atol, double rtol, gsl_integration_method method,
                   double &result, double &error)
{
    switch (method) {
    case romberg: {
        size_t neval;
        const size_t max_workspace_size = 30; // this is GSL's max allowed
        gsl_integration_romberg_workspace * w_romb =
            gsl_integration_romberg_alloc(max_workspace_size);
        gsl_integration_romberg(gsl_func, x_min, x_max, atol, rtol,
                                &result, &neval, w_romb);
        error = GSL_MAX(result*rtol, atol);
        gsl_integration_romberg_free(w_romb);
        return;
        }
    case adaptive_singular: {
        size_t max_partitions = pow(10, 6); // magic number
        gsl_integration_workspace * w_qags =
            gsl_integration_workspace_alloc(max_partitions);
        gsl_integration_qags(gsl_func, x_min, x_max, atol, rtol,
                             max_partitions, w_qags, &result, &error);
        gsl_integration_workspace_free(w_qags);
        return;
        }
    }
    abort();
}
