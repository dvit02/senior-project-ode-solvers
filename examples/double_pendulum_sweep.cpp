//
// Created by Dimitar Vitliyanov on 8.04.26.
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "Integrator.h"
#include "Euler.h"
#include "RK4.h"
#include "DormandPrince.h"
#include "DoublePendulum.h"

// one row of sweep data - everything we want to plot per run
struct SweepRow {
    std::string method;     // "Euler", "AdaptiveRK4", "DP45"
    double tol_or_h;        // rtol for adaptive methods, step size for Euler
    double drift;           // |(E_final - E0) / E0|
    std::size_t rhs_calls;  // total rhs evaluations for this run
    std::size_t steps;      // accepted steps (or total steps for Euler)
};

int main()
{
    // physical parameters - must match the main example
    const double m1 = 1.0;
    const double m2 = 1.0;
    const double L1 = 1.0;
    const double L2 = 1.0;
    const double g  = 9.81;

    DoublePendulum ode(m1, m2, L1, L2, g); // the system
    Integrator integrator; // shared fixed-step / adaptive-RK4 driver

    // time domain - same as main example so results are comparable
    const double t0   = 0.0;
    const double tEnd = 20.0;

    // initial condition - nonzero E0, chaotic regime
    State y0 = {2.0, 2.5, 0.0, 0.0};
    const double E0 = ode.energy(y0); // reference energy for drift

    // accumulator for every run we do
    std::vector<SweepRow> rows;

    // safe relative drift helper - same guard as the main example
    auto rel_drift = [E0](double E) {
        const double denom = std::abs(E0);
        return denom > 1e-12 ? std::abs((E - E0) / E0) : std::abs(E - E0);
    };

    // --- Sweep 1: Euler over decreasing step size ---
    // h values go from loose to tight so we can see convergence
    const std::vector<double> euler_steps = {1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5};

    Euler euler_stepper;

    for (double h : euler_steps)
    {
        Solution sol =
                integrator.integrate(ode, euler_stepper, t0, y0, tEnd, h); // run one Euler pass

        const double E_final = ode.energy(sol.y.back());
        const std::size_t n_steps = sol.t.size() - 1; // every step is kept

        SweepRow row;
        row.method    = "Euler";
        row.tol_or_h  = h;
        row.drift     = rel_drift(E_final);
        row.rhs_calls = n_steps; // Euler: 1 rhs per step
        row.steps     = n_steps;
        rows.push_back(row);

        std::cout << "Euler h=" << h
                  << "  drift=" << row.drift
                  << "  steps=" << n_steps << "\n";
    }

    // --- Sweep 2: Adaptive RK4 over decreasing tolerance ---
    const std::vector<double> tolerances = {1e-4, 1e-6, 1e-8, 1e-10, 1e-15};

    RK4 rk4_stepper;

    for (double rtol : tolerances)
    {
        const double atol  = rtol * 1e-2; // atol two decades tighter - standard ratio
        const double h0    = 0.01;
        const double h_min = 1e-9;
        const double h_max = 0.05;

        AdaptiveResult result =
                integrator.integrateAdaptiveRK4(
                        ode, rk4_stepper,
                        t0, y0, tEnd,
                        h0, rtol, atol,
                        h_min, h_max
                ); // run one adaptive RK4 pass

        const double E_final = ode.energy(result.solution.y.back());
        // step-doubling cost: 3 rk4 calls x 4 evaluations = 12 per trial
        const std::size_t rhs = (result.accepted_steps + result.rejected_steps) * 12;

        SweepRow row;
        row.method    = "AdaptiveRK4";
        row.tol_or_h  = rtol;
        row.drift     = rel_drift(E_final);
        row.rhs_calls = rhs;
        row.steps     = result.accepted_steps;
        rows.push_back(row);

        std::cout << "RK4 rtol=" << rtol
                  << "  drift=" << row.drift
                  << "  rhs=" << rhs << "\n";
    }

    // --- Sweep 3: Dormand-Prince RK45 over decreasing tolerance ---
    DormandPrince dp45;

    for (double rtol : tolerances)
    {
        const double atol  = rtol * 1e-2; // same ratio for fair comparison
        const double h0    = 0.01;
        const double h_min = 1e-9;
        const double h_max = 0.05;

        DP45Result result =
                dp45.integrate(ode, y0, t0, tEnd, h0, rtol, atol, h_min, h_max); // one DP45 pass

        const double E_final = ode.energy(result.solution.y.back());

        SweepRow row;
        row.method    = "DP45";
        row.tol_or_h  = rtol;
        row.drift     = rel_drift(E_final);
        row.rhs_calls = result.rhs_calls; // DP45 tracks this itself
        row.steps     = result.accepted_steps;
        rows.push_back(row);

        std::cout << "DP45 rtol=" << rtol
                  << "  drift=" << row.drift
                  << "  rhs=" << result.rhs_calls << "\n";
    }

    // --- Write everything to one tidy csv ---
    // long format so pandas can pivot easily in the plot script
    std::ofstream out("dpend_sweep.csv");
    out << std::setprecision(12); // enough digits for log-log plots
    out << "method,tol_or_h,drift,rhs_calls,steps\n"; // header

    for (const SweepRow& r : rows)
    {
        out << r.method    << ","
            << r.tol_or_h  << ","
            << r.drift     << ","
            << r.rhs_calls << ","
            << r.steps     << "\n";
    }

    std::cout << "\nWrote " << rows.size() << " rows to dpend_sweep.csv\n";
    return 0;
}