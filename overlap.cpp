

#include "overlap.h"


DetectionCavity::DetectionCavity(double R_init, double L_init,
                                 double seperation_init, double mode_name_init,
                                 VectorFieldInCylinder Ei_mode_init,
                                 CylinderFrequency omega_func_init)
    : R(R_init), L(L_init), seperation(seperation_init),
      mode_name(mode_name_init), Ei_mode(Ei_mode_init),
      omega_func(omega_func_init) {}


class Overlap {
public:
    EffectiveCurrent j_eff;
    DetectionCavity detect_mode;
    double atol, rtol;
    gsl_integration_method method;
    Overlap();
    void operator () (double, double &, double &);
};

Overlap::Overlap()
