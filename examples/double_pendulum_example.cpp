//
// Created by Dimitar Vitliyanov on 8.04.26.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include "Integrator.h"
#include "Euler.h"
#include "RK4.h"
#include "DormandPrince.h"
#include "DoublePendulum.h"

// write a fixed-step or plain Solution with the 4D double pendulum state and energy
static void export_solution(const char* filename,
                            const Solution& sol,
                            const DoublePendulum& ode)
{
    std::ofstream out(filename); // open output file
    out << std::setprecision(12); // enough digits for energy drift analysis in python
    out << "t,theta1,theta2,omega1,omega2,E\n"; // CSV header with energy column

    for (std::size_t i = 0; i < sol.t.size(); ++i) // walk every stored point
    {
        const State& y = sol.y[i];
        out << sol.t[i] << ","
            << y[0] << "," // theta1
            << y[1] << "," // theta2
            << y[2] << "," // omega1
            << y[3] << "," // omega2
            << ode.energy(y) << "\n"; // total mechanical energy
    }
}

// adaptive RK4 wrapper that unpacks the AdaptiveResult
static void export_adaptive(const char* filename,
                            const AdaptiveResult& result,
                            const DoublePendulum& ode)
{
    export_solution(filename, result.solution, ode); // just forward the inner Solution
}

// DP45 wrapper
static void export_dp45(const char* filename,
                        const DP45Result& result,
                        const DoublePendulum& ode)
{
    export_solution(filename, result.solution, ode); // same idea
}

// print one row of the diagnostics table with a label and three numeric columns
static void print_row(const std::string& label,
                      const std::string& euler_val,
                      const std::string& rk4_val,
                      const std::string& dp_val)
{
    const int label_w = 20; // left label column
    const int num_w   = 20; // each numeric column

    std::cout << std::left  << std::setw(label_w) << label
              << std::right << std::setw(num_w)   << euler_val
              << std::right << std::setw(num_w)   << rk4_val
              << std::right << std::setw(num_w)   << dp_val
              << "\n";
}

// format an integer for table display
static std::string fmt_int(std::size_t n)
{
    std::ostringstream os;
    os << n;
    return os.str();
}

// format a floating-point value in scientific notation for tiny drift values
static std::string fmt_sci(double x)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(3) << x;
    return os.str();
}

int main()
{
    // physical parameters - equal masses and rods gives the classic chaotic regime
    const double m1 = 1.0;
    const double m2 = 1.0;
    const double L1 = 1.0;
    const double L2 = 1.0;
    const double g  = 9.81;

    DoublePendulum ode(m1, m2, L1, L2, g); // construct the system
    Integrator integrator; // shared integrator object

    // time domain - 20 seconds is long enough to see chaotic divergence
    const double t0   = 0.0;
    const double tEnd = 20.0;

    // initial condition: both arms lifted but not horizontal - nonzero E0, still chaotic
    State y0 = {2.0, 2.5, 0.0, 0.0};

    // --- Euler (fixed step) ---
    Euler euler_stepper; // explicit Euler stepper
    const double h_euler = 1e-4; // tiny step or Euler will explode energetically

    Solution euler_sol =
            integrator.integrate(ode, euler_stepper, t0, y0, tEnd, h_euler); // run Euler

    export_solution("dpend_euler.csv", euler_sol, ode); // export Euler trajectory

    // --- Adaptive RK4 (step-doubling) ---
    RK4 rk4_stepper; // classical RK4 stepper used inside the adaptive wrapper

    const double h0    = 0.01;  // initial adaptive step size
    const double rtol  = 1e-8;  // tight relative tolerance
    const double atol  = 1e-10; // tight absolute tolerance
    const double h_min = 1e-7;
    const double h_max = 0.05;

    AdaptiveResult rk4_result =
            integrator.integrateAdaptiveRK4(
                    ode, rk4_stepper,
                    t0, y0, tEnd,
                    h0, rtol, atol,
                    h_min, h_max
            ); // run adaptive RK4

    export_adaptive("dpend_rk4.csv", rk4_result, ode); // export adaptive RK4 trajectory

    // --- Dormand-Prince RK45 ---
    DormandPrince dp45; // embedded RK45 integrator

    DP45Result dp_result =
            dp45.integrate(ode, y0, t0, tEnd, h0, rtol, atol, h_min, h_max); // run DP45

    export_dp45("dpend_dp45.csv", dp_result, ode); // export DP45 trajectory

    // --- Diagnostics: side by side ---
    const double E0 = ode.energy(y0); // reference energy for drift comparison

    // safe relative drift helper - falls back to absolute drift when E0 is tiny
    auto rel_drift = [E0](double E) {
        const double denom = std::abs(E0);
        return denom > 1e-12 ? std::abs((E - E0) / E0) : std::abs(E - E0);
    };

    // header block - use scientific for tolerances so they actually show up
    std::cout << "\n";
    std::cout << "Double pendulum  "
              << "m1=" << m1 << "  m2=" << m2
              << "  L1=" << L1 << "  L2=" << L2
              << "  g="  << g  << "\n";
    std::cout << "t in [" << t0 << ", " << tEnd << "]  "
              << std::scientific << std::setprecision(1)
              << "rtol=" << rtol << "  atol=" << atol << "\n";
    std::cout << std::fixed << std::setprecision(6)
              << "Initial energy E0 = " << E0 << "\n\n";

    // table header
    const int total_width = 20 + 20 * 3; // label + three numeric columns
    std::cout << std::string(total_width, '=') << "\n";
    print_row("Quantity", "Euler (fixed)", "Adaptive RK4", "Dormand-Prince");
    std::cout << std::string(total_width, '-') << "\n";

    // steps taken row
    print_row("Steps taken",
              fmt_int(euler_sol.t.size() - 1),           // every Euler step kept
              fmt_int(rk4_result.accepted_steps),        // only accepted RK4 steps
              fmt_int(dp_result.accepted_steps));        // only accepted DP45 steps

    // rejected steps row - Euler has no rejection mechanism
    print_row("Rejected steps",
              "-",
              fmt_int(rk4_result.rejected_steps),
              fmt_int(dp_result.rejected_steps));

    // rhs calls row
    print_row("RHS calls",
              fmt_int(euler_sol.t.size() - 1),           // Euler: 1 rhs per step
              fmt_int((rk4_result.accepted_steps + rk4_result.rejected_steps) * 12),
            // step-doubling RK4: 3 sub-steps x 4 evals = 12 per trial
              fmt_int(dp_result.rhs_calls));

    // final energy drift - tells us how well each method conserves energy
    const double E_euler = ode.energy(euler_sol.y.back());
    const double E_rk4   = ode.energy(rk4_result.solution.y.back());
    const double E_dp    = ode.energy(dp_result.solution.y.back());

    print_row("|dE/E0| final",
              fmt_sci(rel_drift(E_euler)),
              fmt_sci(rel_drift(E_rk4)),
              fmt_sci(rel_drift(E_dp)));

    std::cout << std::string(total_width, '=') << "\n\n";

    return 0;
}