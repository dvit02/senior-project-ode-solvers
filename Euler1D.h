//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_EULER1D_H
#define SENIOR_PROJECT_EULER1D_H
#include "Stepper1D.h"

class Euler1D : public Stepper1D {
public:
    double step(const ODE1D& ode, double t, double y, double h) const override;
};



#endif //SENIOR_PROJECT_EULER1D_H
