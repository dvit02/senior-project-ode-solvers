//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_INTEGRATOR_H
#define SENIOR_PROJECT_INTEGRATOR_H
#include "Solution.h"
#include "Stepper.h"
#include "RK4.h"


struct AdaptiveResult {
    Solution solution;

    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;

    double h_min_used = 0.0;
    double h_max_used = 0.0;
};

class Integrator {
public:
    Solution integrate(const ODE& ode,
                       const Stepper& stepper,
                       double t0,
                       const State& y0,
                       double tEnd,
                       double h) const
    {
        Solution sol;
        sol.t.push_back(t0);
        sol.y.push_back(y0);

        double t = t0;
        State y = y0;
        State y_next;

        while (t < tEnd) {
            double h_step = h;
            if (t + h_step > tEnd) h_step = tEnd - t; // no <algorithm>

            stepper.step(ode, t, y, h_step, y_next);
            t += h_step;
            y = y_next;

            sol.t.push_back(t);
            sol.y.push_back(y);
        }

        return sol;
    }
/**
* Adaptive RK4 integrator using step-doubling error control.
*
* @param ode      ODE system
* @param stepper  RK4 stepper (used internally)
* @param t0       Initial time
* @param y0       Initial state
* @param tEnd     Final time
* @param h0       Initial step size
* @param rtol     Relative tolerance
* @param atol     Absolute tolerance
* @param h_min    Minimum allowed step size
* @param h_max    Maximum allowed step size
*
* @return AdaptiveResult containing solution and diagnostics
     **/


AdaptiveResult integrateAdaptiveRK4(
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
) const;
};
#endif //SENIOR_PROJECT_INTEGRATOR_H
