//
// Created by Dimitar Vitliyanov on 31.03.26.
//
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Euler.h"
#include "ExponentialDecayODE.h"
#include "Integrator.h"
#include "RK4.h"

namespace { // other .cpp files cannot see or accidentally call these functions; modern equivalent of static

    double exact_decay(double t, double y0, double k) // computes the exact solution y(t) = y0 * e^(-kt)
    {
        return y0 * std::exp(-k * t);
    }

    double final_error(const Solution& sol, double y0, double k)
    {
        return std::fabs(sol.y.back()[0] - exact_decay(sol.t.back(), y0, k));
    } // takes numerical solution and exact solution => computes the absolute error at the final time point

    void require(bool condition, const std::string& message)
    {
        if (!condition)
            throw std::runtime_error(message);
    } // if condition is false the program throws a readable message telling you exactly what check failed

    void test_fixed_integrator_hits_terminal_time_and_stores_consistent_grid()
    {
        // parameters
        const double k    = 1.0; // decay rate: y(t) = e^(-t) is the exact solution for this initial condition
        const double t0   = 0.0;
        const double t_end = 1.0;
        const double h    = 0.3; // fixed step size — grid will be {0.0, 0.3, 0.6, 0.9, 1.0}
        const State y0    = {1.0};

        ExponentialDecayODE ode(k);
        Integrator integrator;
        RK4 rk4;

        const Solution sol = integrator.integrate(ode, rk4, t0, y0, t_end, h);
        // runs the full fixed-step integration and stores all time points and states

        // every time point must have a corresponding state and vice versa
        require(sol.t.size() == sol.y.size(),
                "Integrator must store the same number of time and state entries.");

        // h=0.3 on [0,1] gives steps 0.0->0.3->0.6->0.9->1.0 which is exactly 5 points
        require(sol.t.size() == 5,
                "Expected grid points {0.0, 0.3, 0.6, 0.9, 1.0} for h = 0.3 on [0, 1].");

        // the initial condition must be stored before the loop begins
        require(sol.t.front() == t0,
                "First stored time must equal t0.");

        // uses fabs instead of == because the final step is adjusted to hit t_end exactly
        // and floating point arithmetic can introduce tiny rounding errors
        require(std::fabs(sol.t.back() - t_end) < 1e-12,
                "Last stored time must equal t_end after final-step adjustment.");

        for (std::size_t i = 1; i < sol.t.size(); ++i) { // compares each point against the one before it
            require(sol.t[i] > sol.t[i - 1],
                    "Time grid must be strictly increasing."); // t[1] > t[0], t[2] > t[1], etc.
            require(sol.y[i].size() == 1,
                    "Scalar ODE solution should store one component per state."); // this is a 1D ODE
        }
    }

    void test_rk4_is_more_accurate_than_euler_on_the_same_problem()
    {
        // parameters — same step size and interval for both solvers so the comparison is fair
        const double k    = 2.0;
        const double t0   = 0.0;
        const double t_end = 2.0;
        const double h    = 0.2;
        const State y0    = {1.0}; // y(t) = e^(-2t) is the exact solution

        ExponentialDecayODE ode(k);
        Integrator integrator;
        Euler euler; // first-order method
        RK4 rk4;    // fourth-order method

        const Solution euler_sol = integrator.integrate(ode, euler, t0, y0, t_end, h); // Euler integration
        const Solution rk4_sol   = integrator.integrate(ode, rk4,   t0, y0, t_end, h); // RK4 integration

        const double euler_error = final_error(euler_sol, y0[0], k); // Euler absolute error at t_end
        const double rk4_error   = final_error(rk4_sol,   y0[0], k); // RK4 absolute error at t_end

        // RK4 is 4th order and Euler is 1st order so RK4 must always win at the same step size
        require(rk4_error < euler_error,
                "RK4 should be more accurate than Euler on exponential decay for the same step size.");

        // sanity check: even RK4 should not be wildly inaccurate on a smooth scalar ODE
        require(rk4_error < 1e-4,
                "RK4 final-time error is unexpectedly large on a smooth scalar ODE.");
    }

    void test_observed_convergence_rates_follow_theory()
    {
        // convergence rate test: halve h and measure how much the error shrinks
        // theory says Euler error shrinks by ~2x (first order) and RK4 by ~16x (fourth order)
        const double k    = 1.0;
        const double t0   = 0.0;
        const double t_end = 1.0;
        const State y0    = {1.0}; // y(t) = e^(-t) is the exact solution

        ExponentialDecayODE ode(k);
        Integrator integrator;
        Euler euler;
        RK4 rk4;

        const Solution euler_h  = integrator.integrate(ode, euler, t0, y0, t_end, 0.2); // Euler with h
        const Solution euler_h2 = integrator.integrate(ode, euler, t0, y0, t_end, 0.1); // Euler with h/2
        const Solution rk4_h    = integrator.integrate(ode, rk4,   t0, y0, t_end, 0.2); // RK4 with h
        const Solution rk4_h2   = integrator.integrate(ode, rk4,   t0, y0, t_end, 0.1); // RK4 with h/2

        // ratio = error(h) / error(h/2) — how much did the error shrink when we halved the step
        const double euler_ratio = final_error(euler_h,  y0[0], k) / final_error(euler_h2, y0[0], k);
        const double rk4_ratio   = final_error(rk4_h,    y0[0], k) / final_error(rk4_h2,   y0[0], k);

        // first-order method: halving h should halve the error, so ratio ≈ 2
        // window [1.7, 2.3] gives tolerance for the asymptotic regime not being perfectly reached
        require(euler_ratio > 1.7 && euler_ratio < 2.3,
                "Euler should show approximately first-order convergence when h is halved.");

        // fourth-order method: halving h shrinks error by 2^4 = 16, so ratio >> 10
        require(rk4_ratio > 10.0,
                "RK4 should show strong high-order convergence when h is halved.");
    }

} // namespace

int main()
{
    try {
        test_fixed_integrator_hits_terminal_time_and_stores_consistent_grid();
        test_rk4_is_more_accurate_than_euler_on_the_same_problem();
        test_observed_convergence_rates_follow_theory();
    } catch (const std::exception& ex) {
        std::cerr << "Test failure: " << ex.what() << '\n'; // print which test failed and why
        return 1; // non-zero exit code signals failure to the build system
    }

    std::cout << "All fixed-stepper unit tests passed.\n";
    return 0;
}
