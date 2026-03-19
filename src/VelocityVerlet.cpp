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

    // Compute acceleration at current state: a_n = f(t, y)[1]
    State dydt(n);
    ode.rhs(t, y, dydt);
    const double a_n = dydt[1]; // acceleration is second component of RHS

    // Velocity Verlet position update: x_{n+1} = x_n + h*v_n + 0.5*h^2*a_n
    const double x_next = x + h * v + 0.5 * h * h * a_n;

    // Build predicted next state to evaluate acceleration at t_{n+1}
    State y_pred(n);
    y_pred[0] = x_next;
    y_pred[1] = v; // velocity predictor: use current v (standard Verlet convention)

    // Compute acceleration at predicted next state: a_{n+1}
    State dydt_next(n);
    ode.rhs(t + h, y_pred, dydt_next);
    const double a_next = dydt_next[1]; // acceleration at next step

    // Velocity Verlet velocity update: v_{n+1} = v_n + 0.5*h*(a_n + a_{n+1})
    const double v_next = v + 0.5 * h * (a_n + a_next);

    // Store updated state
    y_next[0] = x_next; // updated position
    y_next[1] = v_next; // updated velocity
}