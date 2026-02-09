//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_INTEGRATOR1D_H
#define SENIOR_PROJECT_INTEGRATOR1D_H
#include "Solution1D.h"
#include "Stepper1D.h"

inline Solution1D integrate(const Stepper1D& stepper,
                            const ODE1D& ode,
                            double t0,
                            double y0,
                            double tEnd,
                            double h)
{
    Solution1D sol;

    double t = t0;
    double y = y0;

    sol.t.push_back(t);
    sol.y.push_back(y);

    while (t + 1e-15 < tEnd) {
        y = stepper.step(ode, t, y, h);
        t += h;
        sol.t.push_back(t);
        sol.y.push_back(y);
    }

    return sol;
}

#endif //SENIOR_PROJECT_INTEGRATOR1D_H
