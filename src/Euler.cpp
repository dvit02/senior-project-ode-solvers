//
// Created by Dimitar Vitliyanov on 10.02.26.
//
#include "Euler.h"

void Euler::step(const ODE& ode,
                 double t,
                 const State& y,
                 double h,
                 State& y_next) const
{
    const std::size_t d = ode.dim();

    // Ensure correct size
    if (y_next.size() != d)
        y_next.resize(d);

    State dydt(d);
    ode.rhs(t, y, dydt);

    for (std::size_t i = 0; i < d; ++i) {
        y_next[i] = y[i] + h * dydt[i];
    }
}