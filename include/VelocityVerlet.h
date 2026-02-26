//
// Created by Dimitar Vitliyanov on 26.02.26.
//

#ifndef SENIOR_PROJECT_VELOCITYVERLET_H
#define SENIOR_PROJECT_VELOCITYVERLET_H
#include "Stepper.h"

// Assumes the ODE is in "mechanics form" with y = [x, v]:
// x' = v
// v' = a(t, x, v)
//
// Works for HarmonicOscillator and also VanDerPol (a depends on v).
class VelocityVerlet : public Stepper {
public:
    void step(const ODE& ode,
              double t,
              const State& y,
              double h,
              State& y_next) const override;
};



#endif //SENIOR_PROJECT_VELOCITYVERLET_H
