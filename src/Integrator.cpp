//
// Created by Dimitar Vitliyanov on 25.02.26.
//
#include "Integrator.h"
#include <cmath>
#include <iostream>
AdaptiveResult Integrator::integrateAdaptiveRK4(
        const ODE& ode,
        const RK4& stepper,
        double t0,
        const State& y0,
        double tEnd,
        double h0,
        double rtol,
        double atol,
        double h_min,
        double h_max
) const
{
    AdaptiveResult result;

    double t = t0;
    State y = y0;
    double h = h0;

    const double safety = 0.9;
    const double min_factor = 0.2;
    const double max_factor = 5.0;

    result.solution.t.push_back(t);
    result.solution.y.push_back(y);

    result.h_min_used = h;
    result.h_max_used = h;

    while (t < tEnd)
    {
        if (t + h > tEnd)
            h = tEnd - t;

        // Step-doubling

        State y_big;
        stepper.step(ode, t, y, h, y_big);

        State y_half;
        stepper.step(ode, t, y, 0.5 * h, y_half);

        State y_two;
        stepper.step(ode, t + 0.5 * h, y_half, 0.5 * h, y_two);

        // Error estimation

        double scaled_error = 0.0;

        for (std::size_t i = 0; i < y.size(); ++i)
        {
            double err = std::fabs(y_two[i] - y_big[i]);

            double scale = atol + rtol * (
                    std::fabs(y[i]) > std::fabs(y_two[i])
                    ? std::fabs(y[i])
                    : std::fabs(y_two[i])
            );

            double ratio = err / scale;

            if (ratio > scaled_error)
                scaled_error = ratio;
        }

        //Accept step

        if (scaled_error <= 1.0)
        {
            t += h;
            y = y_two;   // use higher accuracy solution

            result.solution.t.push_back(t);
            result.solution.y.push_back(y);

            result.accepted_steps++;

            if (h < result.h_min_used)
                result.h_min_used = h;

            if (h > result.h_max_used)
                result.h_max_used = h;
        }
        else
        {
            result.rejected_steps++;
        }

        // Step size update

        double factor;

        if (scaled_error == 0.0)
        {
            factor = max_factor;
        }
        else
        {
            factor = safety * std::pow(1.0 / scaled_error, 1.0 / 5.0);
        }

        if (factor < min_factor)
            factor = min_factor;

        if (factor > max_factor)
            factor = max_factor;

        h *= factor;

        if (h < h_min)
            h = h_min;

        if (h > h_max)
            h = h_max;

        if (h <= 0.0)
        {
            std::cerr << "Step size collapsed to zero.\n";
            break;
        }
    }

    return result;
}