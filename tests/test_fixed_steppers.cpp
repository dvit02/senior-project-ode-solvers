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

namespace {

    double exact_decay(double t, double y0, double k)
    {
        return y0 * std::exp(-k * t);
    }

    double final_error(const Solution& sol, double y0, double k)
    {
        return std::fabs(sol.y.back()[0] - exact_decay(sol.t.back(), y0, k));
    }

    void require(bool condition, const std::string& message)
    {
        if (!condition)
            throw std::runtime_error(message);
    }

    void test_fixed_integrator_hits_terminal_time_and_stores_consistent_grid()
    {
        const double k = 1.0;
        const double t0 = 0.0;
        const double t_end = 1.0;
        const double h = 0.3;
        const State y0 = {1.0};

        ExponentialDecayODE ode(k);
        Integrator integrator;
        RK4 rk4;

        const Solution sol = integrator.integrate(ode, rk4, t0, y0, t_end, h);

        require(sol.t.size() == sol.y.size(),
                "Integrator must store the same number of time and state entries.");
        require(sol.t.size() == 5,
                "Expected grid points {0.0, 0.3, 0.6, 0.9, 1.0} for h = 0.3 on [0, 1].");
        require(sol.t.front() == t0,
                "First stored time must equal t0.");
        require(std::fabs(sol.t.back() - t_end) < 1e-12,
                "Last stored time must equal t_end after final-step adjustment.");

        for (std::size_t i = 1; i < sol.t.size(); ++i) {
            require(sol.t[i] > sol.t[i - 1],
                    "Time grid must be strictly increasing.");
            require(sol.y[i].size() == 1,
                    "Scalar ODE solution should store one component per state.");
        }
    }

    void test_rk4_is_more_accurate_than_euler_on_the_same_problem()
    {
        const double k = 2.0;
        const double t0 = 0.0;
        const double t_end = 2.0;
        const double h = 0.2;
        const State y0 = {1.0};

        ExponentialDecayODE ode(k);
        Integrator integrator;
        Euler euler;
        RK4 rk4;

        const Solution euler_sol = integrator.integrate(ode, euler, t0, y0, t_end, h);
        const Solution rk4_sol = integrator.integrate(ode, rk4, t0, y0, t_end, h);

        const double euler_error = final_error(euler_sol, y0[0], k);
        const double rk4_error = final_error(rk4_sol, y0[0], k);

        require(rk4_error < euler_error,
                "RK4 should be more accurate than Euler on exponential decay for the same step size.");
        require(rk4_error < 1e-4,
                "RK4 final-time error is unexpectedly large on a smooth scalar ODE.");
    }

    void test_observed_convergence_rates_follow_theory()
    {
        const double k = 1.0;
        const double t0 = 0.0;
        const double t_end = 1.0;
        const State y0 = {1.0};

        ExponentialDecayODE ode(k);
        Integrator integrator;
        Euler euler;
        RK4 rk4;

        const Solution euler_h   = integrator.integrate(ode, euler, t0, y0, t_end, 0.2);
        const Solution euler_h2  = integrator.integrate(ode, euler, t0, y0, t_end, 0.1);
        const Solution rk4_h     = integrator.integrate(ode, rk4, t0, y0, t_end, 0.2);
        const Solution rk4_h2    = integrator.integrate(ode, rk4, t0, y0, t_end, 0.1);

        const double euler_ratio = final_error(euler_h, y0[0], k) / final_error(euler_h2, y0[0], k);
        const double rk4_ratio   = final_error(rk4_h, y0[0], k) / final_error(rk4_h2, y0[0], k);

        require(euler_ratio > 1.7 && euler_ratio < 2.3,
                "Euler should show approximately first-order convergence when h is halved.");
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
        std::cerr << "Test failure: " << ex.what() << '\n';
        return 1;
    }

    std::cout << "All fixed-stepper unit tests passed.\n";
    return 0;
}