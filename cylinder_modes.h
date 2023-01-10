
#ifndef CYLINDER_MODES_H
#define CYLINDER_MODES_H

#include "cylinder_integration.h"

/* angular frequency of modes */
double angular_frequency_TM010(double, double);

double angular_frequency_TE011(double, double);

/*
Surface current in cylindrical components relative to the (passed)
emitter point. The surface argument indicates on which part of the
cylinder (top, bottom, or sides) the evaluation is occurring. The
two input coordinates (x1, x2) are understood to be cylindrical
(r, phi) for the top and bottom surfaces and (phi, z) for the sides.
*/
enum CylindricalUnitVector {r_hat, phi_hat, z_hat,};

double Ki_cylinder_TM010(double, double, double, double,
                         Surface, CylindricalUnitVector);

double Ki_cylinder_TE011(double, double, double, double,
                         Surface, CylindricalUnitVector);

/*
E field mode in a cylinder cavity, evaluated by cylinder coordinates
with the origin at the center of the cavity endcap.
*/
double Ei_cylinder_TM010(double, double, double,
                         double, double, CylindricalUnitVector);

double Ei_cylinder_TE011(double, double, double,
                         double, double, CylindricalUnitVector);


#endif
