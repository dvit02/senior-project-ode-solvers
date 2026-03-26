//
// Created by Dimitar Vitliyanov on 23.03.26.
//

#ifndef SENIOR_PROJECT_DORMANDPRINCE_H
#define SENIOR_PROJECT_DORMANDPRINCE_H

#include "ODE.h"
#include "State.h"
#include "Solution.h"

// Dormand-Prince RK45 result container
// Mirrors AdaptiveResult but adds RHS call count for benchmarking
struct DP45Result {

    Solution solution;               // computed numerical solution (time points + states)
    std::size_t accepted_steps = 0;  // number of accepted steps
    std::size_t rejected_steps = 0;  // number of rejected steps
    std::size_t rhs_calls = 0;       // total right-hand-side evaluations (6 per trial)
    double h_min_used = 0.0;         // smallest step size used during integration
    double h_max_used = 0.0;         // largest step size used during integration

};

// Dormand-Prince RK45 adaptive integrator
// Uses an embedded 4th/5th order Runge-Kutta pair (Butcher tableau DOPRI5)
// to estimate local error at no extra cost - 6 RHS evaluations per trial
// regardless of acceptance, compared to 12 for step-doubling RK4.
// This is the algorithm behind MATLAB ode45 and SciPy RK45.
class DormandPrince {

public:
    // Integrate ode from t0 to tEnd with adaptive step size control
    DP45Result integrate(const ODE& ode,
                         const State& y0,     // initial state
                         double t0,           // initial time
                         double tEnd,         // final time
                         double h0,           // initial step size
                         double rtol,         // relative tolerance
                         double atol,         // absolute tolerance
                         double h_min,        // minimum allowed step size
                         double h_max         // maximum allowed step size
    ) const;

private:
    // Perform one Dormand-Prince step and return both 4th and 5th order solutions
    // y4: 4th order solution (used for stepping)
    // y5: 5th order solution (used for error estimate only)
    void step(const ODE& ode,
              double t,
              const State& y,
              double h,
              State& y4,     // 4th order output (advance with this)
              State& y5      // 5th order output (error estimate only)
    ) const;
};


#endif //SENIOR_PROJECT_DORMANDPRINCE_H
