//
// Created by Dimitar Vitliyanov on 26.02.26.
//

#include "VelocityVerlet.h"

void VelocityVerlet::step(const ODE& ode,
                          double t,
                          const State& y,
                          double h,
                          State& y_next) const {
    const std::size_t n = ode.dim();
    y_next.assign(n, 0.0);

    // y[0] = position x, y[1] = velocity v
    const double x = y[0];
    const double v = y[1];
}