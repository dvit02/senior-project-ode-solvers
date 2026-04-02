//
// Created by Dimitar Vitliyanov on 2.04.26.
//
#include <iostream>
#include <iomanip>
#include <chrono>
#include "DormandPrince.h"
#include "Integrator.h"
#include "RK4.h"
#include "VanDerPol.h"

using Clock = std::chrono::high_resolution_clock;
using Ms = std::chrono::duration<double, std::milli>;

// Time one run of adaptive RK4 and return the result
static AdaptiveResult run_rk4(const ODE& ode,
                              const RK4& stepper,
                              Integrator& integrator,
                              const State& y0,
                              double t0, double tEnd,
                              double h0, double rtol, double atol,
                              double h_min, double h_max,
                              double& elapsed_ms)
{
    auto start = Clock::now(); // record start time

    AdaptiveResult result =
            integrator.integrateAdaptiveRK4(ode, stepper,
                                            t0, y0, tEnd,
                                            h0, rtol, atol,
                                            h_min, h_max); // run adaptive RK4

    auto end = Clock::now(); // record end time
    elapsed_ms = Ms(end - start).count(); // compute elapsed milliseconds

    return result;
}

// Time one run of Dormand-Prince and return the result
static DP45Result run_dp45(const ODE& ode,
                           DormandPrince& dp45,
                           const State& y0,
                           double t0, double tEnd,
                           double h0, double rtol, double atol,
                           double h_min, double h_max,
                           double& elapsed_ms)
{
    auto start = Clock::now(); // record start time

    DP45Result result =
            dp45.integrate(ode, y0,
                           t0, tEnd,
                           h0, rtol, atol,
                           h_min, h_max); // run Dormand-Prince RK45

    auto end = Clock::now(); // record end time
    elapsed_ms = Ms(end - start).count(); // compute elapsed milliseconds

    return result;
}

int main()
{
    const double mu = 10.0; // Van der Pol stiffness parameter

    VanDerPol ode(mu); // construct Van der Pol ODE system

    // time domain
    const double t0   = 0.0;
    const double tEnd = 20.0;

    State y0 = {2.0, 0.0}; // initial condition: x(0)=2, x'(0)=0

    // shared tolerance settings — identical for both solvers so comparison is fair
    const double h0    = 0.05;  // initial step size
    const double rtol  = 1e-6;  // relative tolerance
    const double atol  = 1e-9;  // absolute tolerance
    const double h_min = 1e-6;  // minimum allowed step size
    const double h_max = 0.1;   // maximum allowed step size

    const int NUM_TRIALS = 5; // number of timing trials (average reduces noise)

   // Adaptive RK4 (step-doubling)
    RK4 rk4_stepper;  // RK4 stepper used inside the adaptive loop
    Integrator integrator;   // integrator that wraps the adaptive step-doubling logic

    AdaptiveResult rk4_result;
    double rk4_total_ms = 0.0;

    for (int i = 0; i < NUM_TRIALS; ++i)
    {
        double elapsed = 0.0;

        rk4_result = run_rk4(ode, rk4_stepper, integrator,
                             y0, t0, tEnd,
                             h0, rtol, atol, h_min, h_max,
                             elapsed); // time one full integration

        rk4_total_ms += elapsed; // accumulate total time
    }

    double rk4_avg_ms = rk4_total_ms / NUM_TRIALS; // average wall time

    // Step-doubling costs 3 RK4 calls per trial step (1 full + 2 half)
    // Each RK4 call makes 4 RHS evaluations  ->  12 RHS per trial step total.
    // This applies to both accepted and rejected steps (no FSAL savings).
    std::size_t rk4_trial_steps = rk4_result.accepted_steps + rk4_result.rejected_steps;
    std::size_t rk4_rhs_calls   = rk4_trial_steps * 12; // estimated RHS evaluations

    //Dormand-Prince RK45
    DormandPrince dp45; // DP45 integrator (embedded 4th/5th order pair)

    DP45Result dp45_result;
    double dp45_total_ms = 0.0;

    for (int i = 0; i < NUM_TRIALS; ++i)
    {
        double elapsed = 0.0;

        dp45_result = run_dp45(ode, dp45,
                               y0, t0, tEnd,
                               h0, rtol, atol, h_min, h_max,
                               elapsed); // time one full integration

        dp45_total_ms += elapsed; // accumulate total time
    }

    double dp45_avg_ms = dp45_total_ms / NUM_TRIALS; // average wall time

    //  Print results
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Van der Pol oscillator benchmark (mu = " << mu << ")\n";
    std::cout << "t in [" << t0 << ", " << tEnd << "]"
              << "  rtol=" << rtol << "  atol=" << atol << "\n";
    std::cout << "Averaged over " << NUM_TRIALS << " trials.\n\n";

    const int w = 26; // column width for table alignment

    std::cout << std::left
              << std::setw(w) << "Quantity"
              << std::setw(w) << "Adaptive RK4"
              << std::setw(w) << "Dormand-Prince RK45"
              << "\n"
              << std::string(w * 3, '-') << "\n";

    std::cout << std::setw(w) << "Wall time (ms)"
              << std::setw(w) << rk4_avg_ms        // average RK4 wall time
              << std::setw(w) << dp45_avg_ms        // average DP45 wall time
              << "\n";

    std::cout << std::setw(w) << "Accepted steps"
              << std::setw(w) << rk4_result.accepted_steps   // steps RK4 accepted
              << std::setw(w) << dp45_result.accepted_steps  // steps DP45 accepted
              << "\n";

    std::cout << std::setw(w) << "Rejected steps"
              << std::setw(w) << rk4_result.rejected_steps   // steps RK4 rejected
              << std::setw(w) << dp45_result.rejected_steps  // steps DP45 rejected
              << "\n";

    std::cout << std::setw(w) << "RHS calls"
              << std::setw(w) << rk4_rhs_calls
              // step-doubling: 12 RHS per trial step (3 RK4 calls x 4 evaluations each)
              << std::setw(w) << dp45_result.rhs_calls
              // DP45 tracks this exactly; includes FSAL reuse across accepted steps
              << "\n";

    std::cout << std::setw(w) << "h_min used"
              << std::setw(w) << rk4_result.h_min_used   // smallest step RK4 needed
              << std::setw(w) << dp45_result.h_min_used  // smallest step DP45 needed
              << "\n";

    std::cout << std::setw(w) << "h_max used"
              << std::setw(w) << rk4_result.h_max_used   // largest step RK4 took
              << std::setw(w) << dp45_result.h_max_used  // largest step DP45 took
              << "\n"
              << std::string(w * 3, '-') << "\n\n";

    // Speed-up ratio: how many times faster is DP45 over adaptive RK4
    double speedup   = rk4_avg_ms / dp45_avg_ms;
    double rhs_ratio = static_cast<double>(rk4_rhs_calls)
                       / static_cast<double>(dp45_result.rhs_calls);

    std::cout << "DP45 speed-up over adaptive RK4: " << speedup   << "x (wall time)\n";
    std::cout << "RK4 uses " << rhs_ratio << "x more RHS calls than DP45\n";

    return 0;
}