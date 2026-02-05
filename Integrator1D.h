//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_INTEGRATOR1D_H
#define SENIOR_PROJECT_INTEGRATOR1D_H
#include "Solution1D.h"


template <typename Solver1D, typename ODEType>
Solution1D integrate(const Solver1D& solver,
                     const ODEType& ode,
                     double t0,
                     double y0,
                     double tEnd,
                     double h)
{
    Solution1D sol;

    double t = t0;
    double y = y0;

    // push initial ONCE
    sol.t.push_back(t);
    sol.y.push_back(y);

    while (t + 1e-15 < tEnd) {
        // advance one step
        y = solver.step(ode, t, y, h);
        t += h;

        // push the new point ONCE
        sol.t.push_back(t);
        sol.y.push_back(y);
    }

    return sol;
}

#endif //SENIOR_PROJECT_INTEGRATOR1D_H
