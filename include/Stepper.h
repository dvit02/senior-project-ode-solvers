//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_STEPPER_H
#define SENIOR_PROJECT_STEPPER_H
#include "ODE.h"

class Stepper {
public:
    virtual ~Stepper() = default;


    virtual void step(const ODE& ode, // Compute one step: y_next ≈ y(t+h)
                      double t,
                      const State& y,
                      double h,
                      State& y_next) const = 0;
};

#endif //SENIOR_PROJECT_STEPPER_H
