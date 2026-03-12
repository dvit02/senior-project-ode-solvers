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
    AdaptiveResult result;      // container for solution and adaptive step statistics

    double t = t0;              // current time
    State y = y0;               // current state
    double h = h0;              // current step size

    const double safety = 0.9;       // safety factor for step size control
    const double min_factor = 0.2;   // minimum step size scaling factor
    const double max_factor = 5.0;   // maximum step size scaling factor

    result.solution.t.push_back(t);      // store initial time
    result.solution.y.push_back(y);      // store initial state

    result.h_min_used = h;      // initialize minimum used step size
    result.h_max_used = h;      // initialize maximum used step size

    while (t < tEnd)
    {
        if (t + h > tEnd)
            h = tEnd - t;       // adjust final step so we stop exactly at tEnd

        // Step-doubling for error estimation

        State y_big;        // one full RK4 step of size h
        stepper.step(ode, t, y, h, y_big);

        State y_half;       // first half-step of size h/2
        stepper.step(ode, t, y, 0.5 * h, y_half);

        State y_two;        // second half-step of size h/2
        stepper.step(ode, t + 0.5 * h, y_half, 0.5 * h, y_two);

        // Error estimation

        double scaled_error = 0.0;      // maximum normalized error over all components

        for (std::size_t i = 0; i < y.size(); ++i)
        {
            double err = std::fabs(y_two[i] - y_big[i]);   // difference between two estimates

            double scale = atol + rtol * (
                    std::fabs(y[i]) > std::fabs(y_two[i])
                    ? std::fabs(y[i])
                    : std::fabs(y_two[i])
            );      // combined absolute and relative tolerance scale

            double ratio = err / scale;      // normalized error in this component

            if (ratio > scaled_error)
                scaled_error = ratio;        // keep the largest componentwise error
        }

        // Accept or reject step

        if (scaled_error <= 1.0)
        {
            t += h;              // accept time advance
            y = y_two;           // use more accurate two-half-step solution

            result.solution.t.push_back(t);      // store accepted time
            result.solution.y.push_back(y);      // store accepted state

            result.accepted_steps++;             // count accepted step

            if (h < result.h_min_used)
                result.h_min_used = h;           // update smallest used step

            if (h > result.h_max_used)
                result.h_max_used = h;           // update largest used step
        }
        else
        {
            result.rejected_steps++;             // count rejected step
        }

        // Step size update

        double factor;       // scaling factor for next step size

        if (scaled_error == 0.0)
        {
            factor = max_factor;     // increase aggressively if error is essentially zero
        }
        else
        {
            factor = safety * std::pow(1.0 / scaled_error, 1.0 / 5.0);
            // RK4-based step update rule using error order
        }

        if (factor < min_factor)
            factor = min_factor;     // do not shrink too aggressively

        if (factor > max_factor)
            factor = max_factor;     // do not grow too aggressively

        h *= factor;         // apply scaling to step size

        if (h < h_min)
            h = h_min;       // enforce minimum step size

        if (h > h_max)
            h = h_max;       // enforce maximum step size

        if (h <= 0.0)
        {
            std::cerr << "Step size collapsed to zero.\n";    // safety check
            break;
        }
    }

    return result;       // return full adaptive integration result
}