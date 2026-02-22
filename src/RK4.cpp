//
// Created by Dimitar Vitliyanov on 19.02.26.
//

#include "RK4.h"
#include "ODE.h"

void RK4::step(const ODE& ode,
               double t,
               double h,
               const std::vector<double>& y,
               std::vector<double>& y_next) const
{
    const std::size_t n = ode.dim();

    // Ensure correct sizes
    y_next.assign(n, 0.0);

    std::vector<double> k1(n), k2(n), k3(n), k4(n);
    std::vector<double> y_tmp(n);

    // k1 = f(t, y)
    ode.rhs(t, y, k1);

    // y_tmp = y + (h/2)*k1
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + 0.5 * h * k1[i];

    // k2 = f(t + h/2, y_tmp)
    ode.rhs(t + 0.5 * h, y_tmp, k2);

    // y_tmp = y + (h/2)*k2
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + 0.5 * h * k2[i];

    // k3 = f(t + h/2, y_tmp)
    ode.rhs(t + 0.5 * h, y_tmp, k3);

    // y_tmp = y + h*k3
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + h * k3[i];

    // k4 = f(t + h, y_tmp)
    ode.rhs(t + h, y_tmp, k4);

    // y_next = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    for (std::size_t i = 0; i < n; ++i)
    {
        y_next[i] = y[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
}

