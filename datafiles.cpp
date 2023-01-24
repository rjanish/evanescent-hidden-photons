#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include "datafiles.h"


double linear_step(double start, double end, size_t N, size_t index)
{
    if (N == 1 && index == 0){
        return start;
    } else if (N > 1 and index < N){
        return start + index*(end - start)/(N - 1);
    } else {
        std::cerr << "invalid sample number" << std::endl;
        abort();
    }
}


double log_step(double start, double end, size_t N, size_t index)
{
    return gsl_sf_exp(
        linear_step(gsl_sf_log(start), gsl_sf_log(end), N, index));
}
