//
// Created by Dimitar Vitliyanov on 23.03.26.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include "DormandPrince.h"
#include "Integrator.h"
#include "RK4.h"
#include "VanDerPol.h"

// Export DP45 solution to CSV
static void export_dp45(const char* filename, const DP45Result& result)
{
    std::ofstream out(filename); // open output file
    out << "t,x,v\n"; // CSV header

    const Solution& sol = result.solution; // reference to solution inside DP45Result

    for (std::size_t i = 0; i < sol.t.size(); ++i) // loop through accepted solution points
    {
        out << sol.t[i] << ","   // write time
            << sol.y[i][0] << ","  // write x
            << sol.y[i][1] << "\n"; // write v
    }
}

// Export adaptive RK4 solution to CSV
static void export_adaptive(const char* filename, const AdaptiveResult& result)
{
    std::ofstream out(filename); // open output file
    out << "t,x,v\n"; // CSV header

    const Solution& sol = result.solution; // reference to solution inside AdaptiveResult

    for (std::size_t i = 0; i < sol.t.size(); ++i) // loop through adaptive solution
    {
        out << sol.t[i] << ","   // write time
            << sol.y[i][0] << ","  // write x
            << sol.y[i][1] << "\n"; // write v
    }
}

int main()
{
    const double mu = 10.0; // Van der Pol parameter

    VanDerPol ode(mu); // construct Van der Pol ODE system

    // time domain
    const double t0   = 0.0;
    const double tEnd = 20.0;

    State y0 = {2.0, 0.0}; // initial condition

    // shared tolerance settings for fair comparison
    const double h0   = 0.05;   // initial step size
    const double rtol = 1e-6;   // relative tolerance
    const double atol = 1e-9;   // absolute tolerance
    const double h_min = 1e-6;  // minimum step size
    const double h_max = 0.1;   // maximum step size

    // --- Dormand-Prince RK45 ---
    DormandPrince dp45; // DP45 integrator

    DP45Result dp_result =
            dp45.integrate(ode, y0, t0, tEnd, h0, rtol, atol, h_min, h_max); // run DP45

    export_dp45("dp45_vdp.csv", dp_result); // export DP45 results

    // --- Adaptive RK4 (step-doubling) for comparison ---
    RK4 rk4_stepper;       // RK4 stepper
    Integrator integrator; // standard integrator

    AdaptiveResult rk4_result =
            integrator.integrateAdaptiveRK4(
                    ode, rk4_stepper,
                    t0, y0, tEnd,
                    h0, rtol, atol,
                    h_min, h_max
            ); // run adaptive RK4 step-doubling

    export_adaptive("rk4_adaptive_vdp.csv", rk4_result); // export adaptive RK4 results

    // --- Diagnostics: side-by-side comparison ---
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Van der Pol oscillator (mu = " << mu << ")\n";
    std::cout << "t in [" << t0 << ", " << tEnd << "]"
              << "  rtol=" << rtol << "  atol=" << atol << "\n\n";

    const int w = 22; // column width

    std::cout << std::left
              << std::setw(w) << "Quantity"
              << std::setw(w) << "Adaptive RK4"
              << std::setw(w) << "Dormand-Prince RK45"
              << "\n"
              << std::string(w * 3, '-') << "\n";

    std::cout << std::setw(w) << "Accepted steps"
              << std::setw(w) << rk4_result.accepted_steps   // RK4 accepted
              << std::setw(w) << dp_result.accepted_steps     // DP45 accepted
              << "\n";

    std::cout << std::setw(w) << "Rejected steps"
              << std::setw(w) << rk4_result.rejected_steps   // RK4 rejected
              << std::setw(w) << dp_result.rejected_steps     // DP45 rejected
              << "\n";

    std::cout << std::setw(w) << "RHS calls"
              << std::setw(w) << (rk4_result.accepted_steps + rk4_result.rejected_steps) * 12
              // step-doubling: 3 RK4 calls x 4 evaluations = 12 per trial
              << std::setw(w) << dp_result.rhs_calls          // DP45 counts directly
              << "\n";

    std::cout << std::setw(w) << "h_min used"
              << std::setw(w) << rk4_result.h_min_used        // RK4 minimum h
              << std::setw(w) << dp_result.h_min_used          // DP45 minimum h
              << "\n";

    std::cout << std::setw(w) << "h_max used"
              << std::setw(w) << rk4_result.h_max_used        // RK4 maximum h
              << std::setw(w) << dp_result.h_max_used          // DP45 maximum h
              << "\n"
              << std::string(w * 3, '-') << "\n";

    return 0;
}