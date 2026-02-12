//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_STEPPER_H
#define SENIOR_PROJECT_STEPPER_H
#include "ODE.h"

// Interface for one-step methods on systems
class Stepper {
public:
    virtual ~Stepper() = default;

    // Compute one step: y_next ≈ y(t+h)
    virtual void step(const ODE& ode,
                      double t,
                      const State& y,
                      double h,
                      State& y_next) const = 0;
};

#endif //SENIOR_PROJECT_STEPPER_H
