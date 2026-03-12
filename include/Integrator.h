//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_INTEGRATOR_H
#define SENIOR_PROJECT_INTEGRATOR_H
#include "Solution.h"
#include "Stepper.h"
#include "RK4.h"

struct AdaptiveResult {

    Solution solution;     // The computed numerical solution (time points + states)
    std::size_t accepted_steps = 0;     // Number of steps accepted by the adaptive algorithm
    std::size_t rejected_steps = 0;     // Number of rejected steps due to large error
    double h_min_used = 0.0;     // Smallest step size used during integration
    double h_max_used = 0.0;     // Largest step size used during integration

};


class Integrator {

public:
    Solution integrate(const ODE& ode,
                       const Stepper& stepper,   // Numerical step method (Euler, RK4, etc.)
                       double t0,                // Initial time
                       const State& y0,          //initial state vector
                       double tEnd,              //final time
                       double h) const           //fixed step size
    {
        Solution sol;        //  solution container


        // Store initial time
        sol.t.push_back(t0);        // Store initial time
        sol.y.push_back(y0);        // -----  initial state

        double t = t0;         // Current time variable
        State y = y0;           // Current state
        State y_next;         // Variable to store next step state


        while (t < tEnd) {         // Loop until final time is reached
            double h_step = h;             // Default step size


            // Adjust step if we would overshoot tEnd
            if (t + h_step > tEnd) h_step = tEnd - t;
            stepper.step(ode, t, y, h_step, y_next);    // Perform one numerical step
            t += h_step;             // Advance time
            y = y_next;            // Update state


            // Store new time, state
            sol.t.push_back(t);
            sol.y.push_back(y);
        }

        // Return the computed solution
        return sol;
    }


    // Adaptive RK4 integrator (error-controlled step size)
    AdaptiveResult integrateAdaptiveRK4(const ODE& ode,
                                        const RK4& stepper,  // RK4 stepper used for integration
                                        double t0,
                                        const State& y0,    // Initial state
                                        double tEnd,
                                        double h0,          // Initial step size
                                        double rtol,        // Relative error tolerance
                                        double atol,        // Absolute error tolerance
                                        double h_min,       // Minimum allowed step size
                                        double h_max        // Maximum allowed step size

    ) const;
};

#endif //SENIOR_PROJECT_INTEGRATOR_H
