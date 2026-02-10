//
// Created by Dimitar Vitliyanov on 9.02.26.
//

#include "../include/RK4_1D.h"

double RK4_1D::step(const ODE1D& ode, double t, double y, double h) const {
    const double k1 = ode.rhs(t, y);
    const double k2 = ode.rhs(t + 0.5*h, y + 0.5*h*k1);
    const double k3 = ode.rhs(t + 0.5*h, y + 0.5*h*k2);
    const double k4 = ode.rhs(t + h,     y + h*k3);

    return y + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}