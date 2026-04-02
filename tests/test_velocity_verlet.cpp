//
// Created by Dimitar Vitliyanov on 31.03.26.
//
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Euler.h"
#include "HarmonicOscillator.h"
#include "Integrator.h"
#include "VelocityVerlet.h"

namespace { // other .cpp files cannot see or accidentally call these functions; modern equivalent of static

    double energy(const State& y, double omega) // computes total mechanical energy: E = 0.5*(v^2 + omega^2 * x^2)
    {
        const double x = y[0]; // position
        const double v = y[1]; // velocity
        return 0.5 * (v * v + omega * omega * x * x); // kinetic + potential energy
    }

    State exact_state(double t, double omega, double x0, double v0)
    {
        const double x = x0 * std::cos(omega * t) + (v0 / omega) * std::sin(omega * t); // exact position
        const double v = -x0 * omega * std::sin(omega * t) + v0 * std::cos(omega * t);  // exact velocity
        return State{x, v}; // exact analytical solution for x'' + omega^2 * x = 0
    }

    double max_energy_deviation(const Solution& sol, double omega)
    {
        const double e0 = energy(sol.y.front(), omega); // initial energy — the ground truth to compare against

        double max_dev = 0.0; // track the worst deviation seen across all time points

        for (const State& y : sol.y)
        {
            const double dev = std::fabs(energy(y, omega) - e0); // how far this point drifted from e0
            if (dev > max_dev)
                max_dev = dev; // keep the largest deviation seen so far
        }

        return max_dev; // returns the worst energy deviation over the entire solution
    }

    void require(bool condition, const std::string& message)
    {
        if (!condition)
            throw std::runtime_error(message);
    } // if condition is false the program throws a readable message telling you exactly what check failed

    void test_velocity_verlet_has_small_long_time_energy_drift()
    {
        // parameters — long integration window to stress-test energy conservation
        const double omega = 1.0; // oscillator frequency
        const double t0    = 0.0;
        const double t_end = 50.0; // 50 time units ~ 8 full oscillation periods
        const double h     = 0.1;  // fixed step size
        const State y0     = {1.0, 0.0}; // x(0)=1, v(0)=0 — starts at maximum displacement

        HarmonicOscillator ode(omega);
        Integrator integrator;
        VelocityVerlet verlet;

        const Solution sol = integrator.integrate(ode, verlet, t0, y0, t_end, h);
        // runs the full fixed-step integration using Velocity Verlet

        const double e0    = energy(sol.y.front(), omega); // initial energy
        const double e_end = energy(sol.y.back(),  omega); // energy at the final time point
        const double drift = std::fabs(e_end - e0);        // total energy drift over the whole run

        // Velocity Verlet is a symplectic integrator — it preserves a modified Hamiltonian
        // so energy should stay bounded and not grow over time the way Euler's does
        require(drift < 2e-3,
                "Velocity Verlet should keep long-time energy drift very small for the harmonic oscillator.");

        // also check the worst deviation at any point mid-integration, not just the endpoint
        require(max_energy_deviation(sol, omega) < 2e-3,
                "Velocity Verlet energy oscillation should remain tightly bounded.");
    }

    void test_velocity_verlet_beats_euler_on_energy_conservation()
    {
        // parameters — same step size and interval for both solvers so the comparison is fair
        const double omega = 1.0;
        const double t0    = 0.0;
        const double t_end = 30.0; // long enough for Euler's energy growth to become clearly visible
        const double h     = 0.1;
        const State y0     = {1.0, 0.0}; // x(0)=1, v(0)=0

        HarmonicOscillator ode(omega);
        Integrator integrator;
        Euler euler;         // first-order, non-symplectic — energy grows over time
        VelocityVerlet verlet; // second-order, symplectic — energy stays bounded

        const Solution euler_sol  = integrator.integrate(ode, euler,  t0, y0, t_end, h); // Euler integration
        const Solution verlet_sol = integrator.integrate(ode, verlet, t0, y0, t_end, h); // Verlet integration

        const double e0           = energy(y0, omega);                                    // exact initial energy
        const double euler_drift  = std::fabs(energy(euler_sol.y.back(),  omega) - e0);  // Euler energy drift at t_end
        const double verlet_drift = std::fabs(energy(verlet_sol.y.back(), omega) - e0);  // Verlet energy drift at t_end

        // Verlet must conserve energy better than Euler — this is the whole point of symplectic integrators
        require(verlet_drift < euler_drift,
                "Velocity Verlet should conserve energy better than Euler on the same oscillator problem.");

        // confirm Euler actually drifts badly — if this fails the test interval is too short to show the effect
        require(euler_drift > 0.1,
                "Euler energy drift should be visibly large in this long-time benchmark.");
    }

    void test_velocity_verlet_tracks_the_exact_solution_over_one_period()
    {
        // integrate exactly one full period T = 2*pi/omega and check position and velocity error
        // after one period the exact solution returns to its initial state: x=1, v=0
        const double omega = 1.0;
        const double t0    = 0.0;
        const double t_end = 2.0 * 3.14159265358979323846; // T = 2*pi for omega=1
        const double h     = 0.01; // small step size so discretisation error stays below 1e-3
        const State y0     = {1.0, 0.0}; // x(0)=1, v(0)=0

        HarmonicOscillator ode(omega);
        Integrator integrator;
        VelocityVerlet verlet;

        const Solution sol   = integrator.integrate(ode, verlet, t0, y0, t_end, h);
        const State exact    = exact_state(sol.t.back(), omega, y0[0], y0[1]); // analytical solution at t_end

        const double position_error = std::fabs(sol.y.back()[0] - exact[0]); // numerical vs exact position
        const double velocity_error = std::fabs(sol.y.back()[1] - exact[1]); // numerical vs exact velocity

        // after one full period both errors should be small — Verlet is second order so h=0.01 is more than enough
        require(position_error < 1e-3,
                "Velocity Verlet position error after one period is too large.");
        require(velocity_error < 1e-3,
                "Velocity Verlet velocity error after one period is too large.");
    }

}

int main()
{
    try {
        test_velocity_verlet_has_small_long_time_energy_drift();
        test_velocity_verlet_beats_euler_on_energy_conservation();
        test_velocity_verlet_tracks_the_exact_solution_over_one_period();
    } catch (const std::exception& ex) {
        std::cerr << "Test failure: " << ex.what() << '\n'; // print which test failed and why
        return 1; // non-zero exit code signals failure to the build system
    }

    std::cout << "All Velocity Verlet unit tests passed.\n";
    return 0;}