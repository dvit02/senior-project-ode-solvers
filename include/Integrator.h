//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_INTEGRATOR_H
#define SENIOR_PROJECT_INTEGRATOR_H
#include "Solution.h"
#include "Stepper.h"

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
};
#endif //SENIOR_PROJECT_INTEGRATOR_H
