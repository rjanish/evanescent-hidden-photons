
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

class CylinderMode {
public:
    double R;
    double L;
    char mode_name[6];
    VectorFieldOnCylinder Ki_emitter;
    CylinderFrequency omega_func;
    double r0, phi0, z0;
    double m;
    Surface surface;
    PropagatorType prop_type;
    CylindricalUnitVector component;
    CylinderMode(double, double, VectorFieldOnCylinder, CylinderFrequency); 
    double Ki_detector(double, double, double);
    double operator () (double, double);
};

class EffectiveCurrent {
public:
    CylinderMode mode; 
    double atol, rtol;
    gsl_integration_method method;
    EffectiveCurrent(double, double, VectorFieldOnCylinder, 
                     CylinderFrequency, double, double,
                     gsl_integration_method);
    void operator () (double, double, double, double, PropagatorType,
                      CylindricalUnitVector, double &, double &);
};


#endif