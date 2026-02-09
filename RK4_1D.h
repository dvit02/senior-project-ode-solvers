//
// Created by Dimitar Vitliyanov on 9.02.26.
//

#ifndef SENIOR_PROJECT_RK4_1D_H
#define SENIOR_PROJECT_RK4_1D_H
#include "Stepper1D.h"

class RK4_1D : public Stepper1D{
public:
    double step(const ODE1D& ode, double t, double y, double h) const override;
};



#endif //SENIOR_PROJECT_RK4_1D_H
