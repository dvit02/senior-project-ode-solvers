//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_EULER_H
#define SENIOR_PROJECT_EULER_H

#include "Stepper.h"

/*Forward Euler method for ODE systems
This declares the step function, which performs one numerical integration step.
It overrides the virtual function from the Stepper base class.
The function takes:
 ode: the differential equation system being solved
t: the current time
y: the current state vector
h: the step size
y_next: the output state after one Euler step
*/
class Euler : public Stepper {
public:
    void step(const ODE& ode,
              double t,
              const State& y,
              double h,
              State& y_next) const override;
};

#endif //SENIOR_PROJECT_EULER_H
