//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#include "../include/Euler1D.h"
double Euler1D::step(const ODE1D& ode, double t, double y, double h) const {
    const double dydt = ode.rhs(t, y);
    return y + h * dydt;
}