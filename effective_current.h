
#ifndef EFFECTIVE_CURRENT_H
#define EFFECTIVE_CURRENT_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>

#include "cylinder_integration.h"
#include "cylinder_modes.h"

#include "effective_current.h"


double wavenumber(double omega, double m);

enum PropagatorType {real, imaginary, evanescent,};

double propagator(double, double, double, double, double,
                  double, double, PropagatorType);

typedef double (*VectorFieldOnCylinder) (double, double, double, double,
                                         Surface, CylindricalUnitVector);

typedef double (*CylinderFrequency) (double, double);

class PropagatedSurfaceCurrent {
public:
    double R;
    double L;
    VectorFieldOnCylinder Ki_emitter;
    CylinderFrequency omega_func;
    double r0, phi0, z0;
    double m;
    Surface surface;
    PropagatorType prop_type;
    CylindricalUnitVector component;
    PropagatedSurfaceCurrent(double, double,
                             VectorFieldOnCylinder,
                             CylinderFrequency);
    double Ki_detector(double, double, double);
    double operator () (double[]);
};

class EffectiveCurrent {
public:
    PropagatedSurfaceCurrent mode;
    double atol, rtol;
    int mineval, maxeval, verbosity;
    EffectiveCurrent(double, double, VectorFieldOnCylinder, 
                     CylinderFrequency, double, double, int, int, int);
    void operator () (double, double, double, double, PropagatorType,
                      CylindricalUnitVector, double &, double &);
};


#endif
