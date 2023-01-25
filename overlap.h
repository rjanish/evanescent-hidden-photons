#ifndef OVERLAP_H
#define OVERLAP_H


#include "cylinder_modes.h"
#include "effective_current.h"


typedef double (*VectorFieldInCylinder) (double, double, double, double,
                                         double, CylindricalUnitVector);


class OverlapIntegrand{
public:
    PropagatedSurfaceCurrent mode;
    double Rd, Ld, seperation;
    VectorFieldInCylinder Edetect;
    double mass;
    PropagatorType re_or_im;
    double atol, rtol;
    int mineval, maxeval, verbosity;
    OverlapIntegrand(double, double, VectorFieldOnCylinder,
                     CylinderFrequency, double, double, double,
                     VectorFieldInCylinder, PropagatorType, double,
                     double, double, int, int, int);
    double operator() (double, double);
};


class Overlap {
public:
    OverlapIntegrand integrand;
    double atol, rtol;
    int mineval, maxeval, verbosity;
    Overlap(double, double, VectorFieldOnCylinder, CylinderFrequency,
            double, double, double, VectorFieldInCylinder,
            double, double, int, int, int);
    double operator() (double);
};

#endif
