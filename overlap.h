#ifndef OVERLAP_H
#define OVERLAP_H


#include "cylinder_modes.h"
#include "effective_current.h"


typedef double (*VectorFieldInCylinder) (double, double, double, double,
                                         double, CylindricalUnitVector);


class OverlapIntegrand{
public:
    EffectiveCurrent j_eff;
    double Rd, Ld, seperation;
    VectorFieldInCylinder Edetect;
    double mass;
    PropagatorType re_or_im;
    double atol, rtol;
    gsl_integration_method method;
    OverlapIntegrand(double, double, VectorFieldOnCylinder,
                     CylinderFrequency, double, double, double,
                     VectorFieldInCylinder, PropagatorType, double,
                     double, double, gsl_integration_method);
    double operator() (double, double);
};


class Overlap {
public:
    OverlapIntegrand integrand;
    double atol, rtol;
    gsl_integration_method method;
    Overlap(double, double, VectorFieldOnCylinder, CylinderFrequency,
            double, double, double, VectorFieldInCylinder,
            double, double, gsl_integration_method);
    double operator() (double);
};

#endif
