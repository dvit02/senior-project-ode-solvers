//
// Created by Dimitar Vitliyanov on 26.02.26.
//
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "Integrator.h"
#include "Euler.h"
#include "RK4.h"
#include "VelocityVerlet.h"
#include "HarmonicOscillator.h"

// Energy for harmonic oscillator: E = 0.5*(v^2 + omega^2 * x^2)
static double energy(const State& y, double omega) {
    const double x = y[0]; // position
    const double v = y[1]; // velocity
    return 0.5 * (v * v + (omega * omega) * (x * x)); // total mechanical energy
}

// Analytical solution: x(t) = x0*cos(wt) + (v0/w)*sin(wt)
static State analytic_state(double t, double omega, double x0, double v0) {
    const double x = x0 * std::cos(omega * t) + (v0 / omega) * std::sin(omega * t); // analytic position
    const double v = -x0 * omega * std::sin(omega * t) + v0 * std::cos(omega * t);  // analytic velocity
    return State{ x, v }; // return analytic state
}

// Export all three methods + analytical to a single CSV
static void export_csv(const char* filename,
                       const Solution& solE,
                       const Solution& solR,
                       const Solution& solV,
                       double omega,
                       double x0,
                       double v0)
{
    std::ofstream out(filename); // open output file
    out << "t,"
           "x_euler,v_euler,E_euler,"
           "x_rk4,v_rk4,E_rk4,"
           "x_verlet,v_verlet,E_verlet,"
           "x_true,v_true,E_true\n"; // CSV header with all methods

    const std::size_t N = std::min({ solE.t.size(), solR.t.size(), solV.t.size() }); // safe minimum

    for (std::size_t i = 0; i < N; ++i) {
        const double t = solE.t[i]; // current time point

        const double xE = solE.y[i][0]; // Euler position
        const double vE = solE.y[i][1]; // Euler velocity
        const double EE = energy(solE.y[i], omega); // Euler energy

        const double xR = solR.y[i][0]; // RK4 position
        const double vR = solR.y[i][1]; // RK4 velocity
        const double ER = energy(solR.y[i], omega); // RK4 energy

        const double xV = solV.y[i][0]; // Verlet position
        const double vV = solV.y[i][1]; // Verlet velocity
        const double EV = energy(solV.y[i], omega); // Verlet energy

        const State yT = analytic_state(t, omega, x0, v0); // exact solution
        const double xT = yT[0]; // analytic position
        const double vT = yT[1]; // analytic velocity
        const double ET = energy(yT, omega); // analytic energy (constant)

        out << t << ","
            << xE << "," << vE << "," << EE << ","
            << xR << "," << vR << "," << ER << ","
            << xV << "," << vV << "," << EV << ","
            << xT << "," << vT << "," << ET << "\n"; // write full row
    }
}

int main() {
    const double omega = 1.0; // oscillator frequency
    HarmonicOscillator ode(omega); // harmonic oscillator ODE system

    Euler euler_stepper;  // forward Euler stepper
    RK4 rk4_stepper;    // classical RK4 stepper
    VelocityVerlet verlet_stepper; // symplectic Velocity Verlet stepper

    Integrator integrator; // integrator controlling time stepping

    const double t0 = 0.0;   // start time
    const double tEnd = 50.0;  // long integration to reveal energy drift
    const double h = 0.1;   // fixed step size for all methods

    State y0 = {1.0, 0.0}; // IC: x(0)=1, v(0)=0

    // Integrate all three methods on the same grid
    const Solution solE = integrator.integrate(ode, euler_stepper, t0, y0, tEnd, h); // Euler
    const Solution solR = integrator.integrate(ode, rk4_stepper, t0, y0, tEnd, h); // RK4
    const Solution solV = integrator.integrate(ode, verlet_stepper, t0, y0, tEnd, h); // Verlet

    // Export combined CSV for Python plotting
    export_csv("verlet_comparison.csv", solE, solR, solV, omega, y0[0], y0[1]);
    std::cout << "Saved CSV: verlet_comparison.csv\n\n";

    // Console output: energy at start, midpoint, end
    const std::size_t N = std::min({solE.t.size(), solR.t.size(), solV.t.size()});
    const std::size_t mid = N / 2; // midpoint index

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\nEnergy conservation comparison"
              << "  (omega=" << omega
              << ", h=" << h
              << ", t in [" << t0 << ", " << tEnd << "])\n\n";

    const int wt = 10;
    const int we = 16;

    std::cout << std::right
              << std::setw(wt) << "t"
              << std::setw(we) << "E_Euler"
              << std::setw(we) << "E_RK4"
              << std::setw(we) << "E_Verlet"
              << std::setw(we) << "E_exact"
              << "\n"
              << std::string(wt + we * 4, '-') << "\n";

    for (std::size_t i : { std::size_t(0), mid, N - 1 }) {
        const double t_i = solE.t[i];
        const State  yT  = analytic_state(t_i, omega, y0[0], y0[1]);

        std::cout << std::right
                  << std::setw(wt) << std::fixed << std::setprecision(2) << t_i
                  << std::setw(we) << std::fixed << std::setprecision(8) << energy(solE.y[i], omega)
                  << std::setw(we) << std::fixed << std::setprecision(8) << energy(solR.y[i], omega)
                  << std::setw(we) << std::fixed << std::setprecision(8) << energy(solV.y[i], omega)
                  << std::setw(we) << std::fixed << std::setprecision(8) << energy(yT, omega)
                  << "\n";
    }

    std::cout << std::string(wt + we * 4, '-') << "\n";
    return 0;
}