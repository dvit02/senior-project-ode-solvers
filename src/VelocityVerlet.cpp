//
// Created by Dimitar Vitliyanov on 26.02.26.
//

#include "VelocityVerlet.h"

void VelocityVerlet::step(const ODE& ode,
                          double t,
                          const State& y,
                          double h,
                          State& y_next) const
{
    const std::size_t n = ode.dim();
    y_next.assign(n, 0.0);

    // y[0] = position x, y[1] = velocity v
    const double x = y[0];
    const double v = y[1];

    // Compute acceleration at current position: a(t)
    State dydt(n);
    ode.rhs(t, y, dydt);
    const double a = dydt[1];

    // Update position: x(t+h) = x + v*h + 0.5*a*h^2
    const double x_next = x + v * h + 0.5 * a * h * h;

    // Compute acceleration at new position: a(t+h)
    State y_mid(n);
    y_mid[0] = x_next;
    y_mid[1] = v;
    State dydt_next(n);
    ode.rhs(t + h, y_mid, dydt_next);
    const double a_next = dydt_next[1];

    // Update velocity: v(t+h) = v + 0.5*(a + a_next)*h
    const double v_next = v + 0.5 * (a + a_next) * h;

    y_next[0] = x_next;
    y_next[1] = v_next;
}