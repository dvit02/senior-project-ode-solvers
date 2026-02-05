//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_EULER1D_H
#define SENIOR_PROJECT_EULER1D_H
#include "ODE1D.h"
class Euler1D {
public:
    // One Euler step: (t, y) -> (t+h, y_next)
    double step(const ODE1D& ode, double t, double y, double h) const;
};



#endif //SENIOR_PROJECT_EULER1D_H
