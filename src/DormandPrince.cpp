//
// Created by Dimitar Vitliyanov on 23.03.26.
//

#include "DormandPrince.h"
#include <cmath>
#include <iostream>

// Dormand-Prince Butcher tableau coefficients (DOPRI5)
// Six stages: c2..c6 are time nodes, a21..a65 are stage weights,
// b4[0..6] are 4th order weights, b5[0..6] are 5th order weights.
// The 5th order solution reuses all six stages - error estimate is free.

// Time nodes
static const double c2 = 1.0 / 5.0;
static const double c3 = 3.0 / 10.0;
static const double c4 = 4.0 / 5.0;
static const double c5 = 8.0 / 9.0;
// c6 = 1.0, c7 = 1.0 (not needed explicitly)

// Stage weights (row i of Butcher tableau)
static const double a21 = 1.0 / 5.0;

static const double a31 = 3.0 / 40.0;
static const double a32 = 9.0 / 40.0;

static const double a41 =  44.0 / 45.0;
static const double a42 = -56.0 / 15.0;
static const double a43 =  32.0 / 9.0;

static const double a51 =  19372.0 / 6561.0;
static const double a52 = -25360.0 / 2187.0;
static const double a53 =  64448.0 / 6561.0;
static const double a54 =   -212.0 /  729.0;

static const double a61 =   9017.0 / 3168.0;
static const double a62 =   -355.0 /   33.0;
static const double a63 =  46732.0 / 5247.0;
static const double a64 =     49.0 /  176.0;
static const double a65 =  -5103.0 / 18656.0;

// 4th order solution weights (advance with these)
static const double b4_1 =   35.0 /  384.0;
// b4_2 = 0
static const double b4_3 =  500.0 / 1113.0;
static const double b4_4 =  125.0 /  192.0;
static const double b4_5 = -2187.0 / 6784.0;
static const double b4_6 =   11.0 /   84.0;
// b4_7 = 0

// 5th order solution weights (error estimate only)
static const double b5_1 =  5179.0 / 57600.0;
// b5_2 = 0
static const double b5_3 =  7571.0 / 16695.0;
static const double b5_4 =   393.0 /   640.0;
static const double b5_5 = -92097.0 / 339200.0;
static const double b5_6 =   187.0 /  2100.0;
static const double b5_7 =     1.0 /    40.0;

void DormandPrince::step(const ODE& ode,
                         double t,
                         const State& y,
                         double h,
                         State& y4,
                         State& y5) const
{
    const std::size_t n = ode.dim();

    // Allocate all six stage derivative vectors
    State k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n);
    State y_tmp(n); // temporary state for stage evaluations

    // k1 = f(t, y)
    ode.rhs(t, y, k1);

    // k2 = f(t + c2*h,  y + h*a21*k1)
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + h * a21 * k1[i];
    ode.rhs(t + c2 * h, y_tmp, k2);

    // k3 = f(t + c3*h,  y + h*(a31*k1 + a32*k2))
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
    ode.rhs(t + c3 * h, y_tmp, k3);

    // k4 = f(t + c4*h,  y + h*(a41*k1 + a42*k2 + a43*k3))
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    ode.rhs(t + c4 * h, y_tmp, k4);

    // k5 = f(t + c5*h,  y + h*(a51*k1 + ... + a54*k4))
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i]
                               + a53 * k3[i] + a54 * k4[i]);
    ode.rhs(t + c5 * h, y_tmp, k5);

    // k6 = f(t + h,  y + h*(a61*k1 + ... + a65*k5))
    for (std::size_t i = 0; i < n; ++i)
        y_tmp[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i]
                               + a64 * k4[i] + a65 * k5[i]);
    ode.rhs(t + h, y_tmp, k6);

    // 4th order solution: y4 = y + h*(b4_1*k1 + b4_3*k3 + b4_4*k4 + b4_5*k5 + b4_6*k6)
    // Note: b4_2 = 0 and b4_7 = 0 so k2 and k7 drop out
    y4.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        y4[i] = y[i] + h * (b4_1 * k1[i] + b4_3 * k3[i]
                            + b4_4 * k4[i] + b4_5 * k5[i] + b4_6 * k6[i]);

    // k7 = f(t + h,  y4)  — FSAL: first same as last
    // This evaluation is reused as k1 of the next accepted step
    ode.rhs(t + h, y4, k7);

    // 5th order solution: y5 = y + h*(b5_1*k1 + b5_3*k3 + ... + b5_7*k7)
    // Note: b5_2 = 0 so k2 drops out
    y5.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        y5[i] = y[i] + h * (b5_1 * k1[i] + b5_3 * k3[i] + b5_4 * k4[i]
                            + b5_5 * k5[i] + b5_6 * k6[i] + b5_7 * k7[i]);
}

DP45Result DormandPrince::integrate(const ODE& ode,
                                    const State& y0,
                                    double t0,
                                    double tEnd,
                                    double h0,
                                    double rtol,
                                    double atol,
                                    double h_min,
                                    double h_max) const
{
    DP45Result result;      // container for solution and diagnostics

    double t = t0;          // current time
    State  y = y0;          // current state
    double h = h0;          // current step size

    const double safety     = 0.9;   // safety factor for step size control
    const double min_factor = 0.2;   // minimum step size scaling factor
    const double max_factor = 5.0;   // maximum step size scaling factor

    result.solution.t.push_back(t);   // store initial time
    result.solution.y.push_back(y);   // store initial state

    result.h_min_used = h;   // initialize minimum used step size
    result.h_max_used = h;   // initialize maximum used step size

    while (t < tEnd)
    {
        if (t + h > tEnd)
            h = tEnd - t;    // adjust final step so we land exactly on tEnd

        State y4, y5;        // 4th and 5th order solutions for this trial
        step(ode, t, y, h, y4, y5);

        result.rhs_calls += 6;   // 6 RHS evaluations per trial (+ 1 FSAL reuse)

        // Error estimation: componentwise scaled difference between 4th and 5th order
        double scaled_error = 0.0;   // maximum normalised error over all components

        for (std::size_t i = 0; i < y.size(); ++i)
        {
            double err = std::fabs(y5[i] - y4[i]);   // difference between 4th and 5th order

            double scale = atol + rtol * (
                    std::fabs(y[i]) > std::fabs(y4[i])
                    ? std::fabs(y[i])
                    : std::fabs(y4[i])
            );   // combined absolute and relative tolerance scale

            double ratio = err / scale;   // normalised error in this component

            if (ratio > scaled_error)
                scaled_error = ratio;    // keep the largest componentwise error
        }

        // Accept or reject step

        if (scaled_error <= 1.0)
        {
            t += h;     // accept time advance
            y = y4;     // advance with 4th order solution

            result.solution.t.push_back(t);   // store accepted time
            result.solution.y.push_back(y);   // store accepted state

            result.accepted_steps++;          // count accepted step

            if (h < result.h_min_used)
                result.h_min_used = h;        // update smalest used step

            if (h > result.h_max_used)
                result.h_max_used = h;        // update largest used step
        }
        else
        {
            result.rejected_steps++;          // count rejected step
        }

        // Step size update using RK45 optimal controller (exponent 1/5)

        double factor;   // scaling factor for next step size

        if (scaled_error == 0.0)
        {
            factor = max_factor;   // increase aggressively if error is essentially zero
        }
        else
        {
            factor = safety * std::pow(1.0 / scaled_error, 1.0 / 5.0);
            // DP45-based step update rule using embedded error order
        }

        if (factor < min_factor)
            factor = min_factor;   // do not shrink too aggressively

        if (factor > max_factor)
            factor = max_factor;   // do not grow too aggressively

        h *= factor;   // apply scaling to step size

        if (h < h_min)
            h = h_min;   // enforce minimum step size

        if (h > h_max)
            h = h_max;   // enforce maximum step size

        if (h <= 0.0)
        {
            std::cerr << "Step size collapsed to zero.\n";   // safety check
            break;
        }
    }

    return result;
}
