//
// Created by Dimitar Vitliyanov on 9.02.26.
//

#ifndef SENIOR_PROJECT_STEPPER1D_H
#define SENIOR_PROJECT_STEPPER1D_H
#include "ODE1D.h"

class Stepper1D {
public:
    virtual ~Stepper1D() = default;

    virtual double step(const ODE1D& ode,
                        double t,
                        double y,
                        double h) const = 0;
};
#endif //SENIOR_PROJECT_STEPPER1D_H
