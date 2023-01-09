#ifndef OVERLAP_H
#define OVERLAP_H


#include "cylinder_modes.h"
#include "effective_current.h"


typedef double (*VectorFieldInCylinder) (double, double, double, double,
                                         double, CylindricalUnitVector));

class DetectionCavity {
public:
    double R;
    double L;
    double seperation;
    char mode_name[6];
    VectorFieldInCylinder Ei_mode;
    CylinderFrequency omega_func;
    DetectionCavity()
};


class Overlap {
public:
    EffectiveCurrent j_eff;
    DetectionCavity detect_mode;
    double atol, rtol;
    gsl_integration_method method;
    Overlap();
    void operator () (double, double &, double &);
};

#endif
