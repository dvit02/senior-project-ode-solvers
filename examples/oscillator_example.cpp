//
// Created by Dimitar Vitliyanov on 10.02.26.
//
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "Euler.h"
#include "RK4.h"
#include "Integrator.h"
#include "Solution.h"
#include "State.h"
#include "HarmonicOscillator.h"

// Energy for harmonic oscillator: E = 0.5*(v^2 + omega^2 x^2)
static double energy(const State& y, double omega) {
    const double x = y[0];
    const double v = y[1];
    return 0.5 * (v*v + (omega*omega) * (x*x));
}

// Analytical solution for x'' + omega^2 x = 0 with IC: x(0)=x0, v(0)=v0
static State analytic_state(double t, double omega, double x0, double v0) {
    const double x = x0 * std::cos(omega * t) + (v0 / omega) * std::sin(omega * t);
    const double v = -x0 * omega * std::sin(omega * t) + v0 * std::cos(omega * t);
    return State{ x, v };
}

// One combined CSV: Euler + RK4 + Analytical ("true") on the same time grid
static void export_csv_combined(const char* filename,
                                const Solution& solE,
                                const Solution& solR,
                                double omega,
                                double x0,
                                double v0)
{
    std::ofstream out(filename);
    out << "t,"
           "x_euler,v_euler,E_euler,"
           "x_rk4,v_rk4,E_rk4,"
           "x_true,v_true,E_true\n";

    const std::size_t N = std::min(solE.t.size(), solR.t.size());

    for (std::size_t i = 0; i < N; ++i) {
        const double t  = solE.t[i];

        const double xE = solE.y[i][0];
        const double vE = solE.y[i][1];
        const double EE = energy(solE.y[i], omega);

        const double xR = solR.y[i][0];
        const double vR = solR.y[i][1];
        const double ER = energy(solR.y[i], omega);

        const State yT  = analytic_state(t, omega, x0, v0);
        const double xT = yT[0];
        const double vT = yT[1];
        const double ET = energy(yT, omega);

        out << t << ","
            << xE << "," << vE << "," << EE << ","
            << xR << "," << vR << "," << ER << ","
            << xT << "," << vT << "," << ET << "\n";
    }
}

int main() {
    const double omega = 1.0;
    HarmonicOscillator ode(omega);

    Euler euler_stepper;
    RK4   rk4_stepper;

    Integrator integrator;

    const double t0   = 0.0;
    const double tEnd = 10.0;
    const double h    = 0.01;

    State y0 = {1.0, 0.0}; // x(0)=1, v(0)=0

    // Print interval (DISPLAY only)
    const double print_every_dt = 2.0;
    const std::size_t print_every_steps =
            static_cast<std::size_t>(std::max<long long>(1LL, std::llround(print_every_dt / h)));

    // Integrate BOTH methods on the SAME grid
    const Solution solE = integrator.integrate(ode, euler_stepper, t0, y0, tEnd, h);
    const Solution solR = integrator.integrate(ode, rk4_stepper,   t0, y0, tEnd, h);

    // Safety checks
    if (solE.t.empty() || solE.y.empty()) {
        std::cout << "Euler solution is empty (no steps produced).\n";
        return 0;
    }
    if (solR.t.empty() || solR.y.empty()) {
        std::cout << "RK4 solution is empty (no steps produced).\n";
        return 0;
    }

    // Ensure same size
    const std::size_t N = std::min(solE.t.size(), solR.t.size());

    // Export ONE combined CSV for Python plotting (Euler + RK4 + Analytical)
    export_csv_combined("oscillator_validation.csv", solE, solR, omega, y0[0], y0[1]);
    std::cout << "Saved CSV: oscillator_validation.csv\n\n";

    // Output header
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Harmonic Oscillator: Explicit Euler vs RK4\n";
    std::cout << "System (ODE): x'(t)=v(t), v'(t)=-(omega^2)x(t)\n";
    std::cout << "omega = " << omega << "\n";
    std::cout << "IC: x(" << t0 << ")=" << y0[0] << ", v(" << t0 << ")=" << y0[1] << "\n";
    std::cout << "h = " << h << ", t in [" << t0 << ", " << tEnd << "]\n\n";

    std::cout << "Energy: E(t)=0.5*(v^2 + omega^2 x^2). Exact model conserves energy.\n";
    std::cout << "Printing every " << print_every_steps << " steps (≈ "
              << (static_cast<double>(print_every_steps) * h) << " time units)\n\n";

    // Layout: n | t | Euler: x v E | RK4: x v E | |E_E - E_R|
    std::cout << std::left
              << std::setw(8)  << "n"
              << std::setw(12) << "t"
              << std::setw(16) << "x_Euler"
              << std::setw(16) << "v_Euler"
              << std::setw(18) << "E_Euler"
              << "  "
              << std::setw(16) << "x_RK4"
              << std::setw(16) << "v_RK4"
              << std::setw(18) << "E_RK4"
              << "  "
              << std::setw(16) << "|E_E-E_R|"
              << "\n";

    std::cout << std::string(8+12+16+16+18 +2 +16+16+18 +2 +16, '-') << "\n";

    auto print_row = [&](std::size_t i) {
        const double t = solE.t[i];
        const std::size_t n = static_cast<std::size_t>(std::llround((t - t0) / h));

        const double xE = solE.y[i][0];
        const double vE = solE.y[i][1];
        const double EE = energy(solE.y[i], omega);

        const double xR = solR.y[i][0];
        const double vR = solR.y[i][1];
        const double ER = energy(solR.y[i], omega);

        std::cout << std::setw(8)  << n
                  << std::setw(12) << t
                  << std::setw(16) << xE
                  << std::setw(16) << vE
                  << std::setw(18) << EE
                  << "  "
                  << std::setw(16) << xR
                  << std::setw(16) << vR
                  << std::setw(18) << ER
                  << "  "
                  << std::setw(16) << std::fabs(EE - ER)
                  << "\n";
    };

    // Print first
    print_row(0);

    // Print every K steps (avoid printing last twice)
    for (std::size_t i = print_every_steps; i < N; i += print_every_steps) {
        if (i + print_every_steps >= N) break;
        print_row(i);
    }

    // Print last
    if (N > 1) print_row(N - 1);

    return 0;
}