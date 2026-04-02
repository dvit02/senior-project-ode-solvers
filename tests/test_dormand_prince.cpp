//
// Created by Dimitar Vitliyanov on 31.03.26.
//
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "DormandPrince.h"
#include "ExponentialDecayODE.h"
#include "Integrator.h"
#include "RK4.h"

namespace { //  other .cpp files cannot see or accidentally call these functions; modern equivelent of static

    double exact_decay(double t, double y0, double k) //  computes the exact solution
    {
        return y0 * std::exp(-k * t);
    }

    double final_error(const Solution& sol, double y0, double k)
    {
        const double t = sol.t.back();
        const double exact = exact_decay(t, y0, k);
        return std::fabs(sol.y.back()[0] - exact);
    } // takes numerical solution and exact solution => computes the error

    void require(bool condition, const std::string& message)
    {
        if (!condition)
            throw std::runtime_error(message);
    }// if program crashes it throws a readeble messege

    void test_dp45_preserves_core_invariants() // unit test begins
    {
        // parameters
        const double k = 1.5;
        const double t0 = 0.0;
        const double t_end = 4.0;
        const double h0 = 0.35;
        const double rtol = 1e-6;
        const double atol = 1e-9;
        const double h_min = 1e-5;
        const double h_max = 0.5;
        const State y0 = {1.0}; // y(t) = e^(-1.5t). is the excat soution for this initial condition

        ExponentialDecayODE ode(k);
        DormandPrince solver;

        const DP45Result result =
                solver.integrate(ode, y0, t0, t_end, h0, rtol, atol, h_min, h_max);
        // This constructs the ODE and the solver, runs the full integration, and stores everything

        require(!result.solution.t.empty(), "DP45 returned an empty time grid.");

        // Every time point must have a corresponding state, and every state must have a corresponding time.
        require(result.solution.t.size() == result.solution.y.size(),
                "Time and state arrays must have the same length.");

        require(result.solution.t.front() == t0,
                "First stored time must equal the initial time.");

        require(std::fabs(result.solution.t.back() - t_end) < 1e-12,
                "Last stored time must land exactly on t_end.");
        // because we first sotre the started stpe
        require(result.accepted_steps + 1 == result.solution.t.size(),
                "Accepted step count must match the number of stored solution points.");

        require(result.rhs_calls == 6 * (result.accepted_steps + result.rejected_steps),
                "Dormand-Prince RHS accounting is inconsistent with 6 calls per trial step.");

        require(result.accepted_steps > 0,
                "Adaptive integration should accept at least one step.");

        require(result.h_min_used >= h_min,
                "Reported minimum used step size is below the configured lower bound.");

        require(result.h_max_used <= h_max,
                "Reported maximum used step size is above the configured upper bound.");

        require(result.h_max_used > result.h_min_used,
                "Adaptive solver did not vary the step size on a nontrivial problem.");

        for (std::size_t i = 1; i < result.solution.t.size(); ++i) { // compares each point against the one before it
            require(result.solution.t[i] > result.solution.t[i - 1],
                    "Accepted time points must be strictly increasing."); // t[1] > t[0]
            require(result.solution.y[i][0] > 0.0,
                    "Exponential decay solution should remain positive for positive initial data.");
        }
    }

    void test_tighter_tolerance_improves_dp45_accuracy()
    {
        const double k = 2.0;
        const double t0 = 0.0;
        const double t_end = 5.0;
        const double h0 = 0.4;
        const double h_min = 1e-6;
        const double h_max = 0.75;
        const State y0 = {1.0};

        ExponentialDecayODE ode(k);
        DormandPrince solver;

        const DP45Result loose =
                solver.integrate(ode, y0, t0, t_end, h0, 1e-3, 1e-6, h_min, h_max);

        const DP45Result tight =
                solver.integrate(ode, y0, t0, t_end, h0, 1e-8, 1e-11, h_min, h_max);

        const double loose_error = final_error(loose.solution, y0[0], k);
        const double tight_error = final_error(tight.solution, y0[0], k);

        require(tight_error < loose_error,
                "Tighter tolerances should reduce the final-time global error.");
        require(tight.accepted_steps >= loose.accepted_steps,
                "Tighter tolerances should not require fewer accepted steps.");
        require(tight.h_min_used <= loose.h_min_used,
                "Tighter tolerances should allow steps at least as small as the loose run.");
    }

    void test_dp45_matches_adaptive_rk4_on_exact_solution()
    {
        const double k = 1.0;
        const double t0 = 0.0;
        const double t_end = 3.0;
        const double h0 = 0.25;
        const double rtol = 1e-6;
        const double atol = 1e-9;
        const double h_min = 1e-6;
        const double h_max = 0.5;
        const State y0 = {1.0};

        ExponentialDecayODE ode(k);
        DormandPrince dp45;
        Integrator integrator;
        RK4 rk4;

        const DP45Result dp_result =
                dp45.integrate(ode, y0, t0, t_end, h0, rtol, atol, h_min, h_max);

        const AdaptiveResult rk4_result =
                integrator.integrateAdaptiveRK4(ode, rk4, t0, y0, t_end, h0, rtol, atol, h_min, h_max);

        const double exact = exact_decay(t_end, y0[0], k);
        const double dp_error = std::fabs(dp_result.solution.y.back()[0] - exact);
        const double rk4_error = std::fabs(rk4_result.solution.y.back()[0] - exact);

        require(dp_error < 1e-5,
                "Dormand-Prince final value is not accurate enough on exponential decay.");
        // DP45's final answer must be within 1e-5 of the true value e^(-1.0 * 3.0).

        require(rk4_error < 1e-5,
                "Adaptive RK4 final value is not accurate enough on exponential decay.");
        //Adaptive final answer must be within 1e-5 of the true value e^(-1.0 * 3.0).

        require(std::fabs(dp_result.solution.y.back()[0] - rk4_result.solution.y.back()[0]) < 1e-5,
                "Dormand-Prince and adaptive RK4 should agree closely on a smooth scalar test problem.");
    }

}

int main()
{
    try {
        test_dp45_preserves_core_invariants();
        test_tighter_tolerance_improves_dp45_accuracy();
        test_dp45_matches_adaptive_rk4_on_exact_solution();
    } catch (const std::exception& ex) {
        std::cerr << "Test failure: " << ex.what() << '\n';
        return 1;
    }

    std::cout << "All Dormand-Prince unit tests passed.\n";
    return 0;
}
