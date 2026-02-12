//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_EULER_H
#define SENIOR_PROJECT_EULER_H

#include "Stepper.h"

// Forward Euler method for ODE systems
class Euler : public Stepper {
public:
    void step(const ODE& ode,
              double t,
              const State& y,
              double h,
              State& y_next) const override;
};

#endif //SENIOR_PROJECT_EULER_H
