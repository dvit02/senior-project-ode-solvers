#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>   // std::max
#include <cstddef>

#include "Euler.h"
#include "Integrator.h"
#include "ODE.h"
#include "Solution.h"
#include "State.h"

// Harmonic oscillator: x' = v,  v' = -x
class HarmonicOscillator : public ODE {
public:
    std::size_t dim() const override { return 2; }

    void rhs(double /*t*/, const State& y, State& dydt) const override {
        dydt.resize(2);
        const double x = y[0];
        const double v = y[1];
        dydt[0] = v;
        dydt[1] = -x;
    }
};

static double energy(const State& y) {
    const double x = y[0];
    const double v = y[1];
    return 0.5 * (x*x + v*v);
}

int main() {
    HarmonicOscillator ode;
    Euler stepper;
    Integrator integrator;

    // ---- Settings ----
    const double t0 = 0.0;
    const double tEnd = 10.0;
    const double h = 0.01;

    State y0 = {1.0, 0.0};   // x(0)=1, v(0)=0

    // Print interval (DISPLAY only)
    const double print_every_dt = 2.0;
    const std::size_t print_every_steps =
            static_cast<std::size_t>(std::max(1.0, std::round(print_every_dt / h)));

    // Integrate
    const Solution sol = integrator.integrate(ode, stepper, t0, y0, tEnd, h);

    // Safety check
    if (sol.t.empty() || sol.y.empty()) {
        std::cout << "Solution is empty (no steps produced).\n";
        return 0;
    }

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Harmonic Oscillator solved with Explicit Euler\n";
    std::cout << "ODE system:\n";
    std::cout << "  x'(t) = v(t)\n";
    std::cout << "  v'(t) = -x(t)\n";
    std::cout << "Initial conditions:\n";
    std::cout << "  x(" << t0 << ") = " << y0[0] << ",  v(" << t0 << ") = " << y0[1] << "\n\n";

    std::cout << "Step size (time increment): h = " << h << "\n";
    std::cout << "Meaning of columns:\n";
    std::cout << "  n        = step index (number of Euler updates performed)\n";
    std::cout << "  t (time) = physical time, computed as  t_n = t0 + n*h\n";
    std::cout << "  x(t)     = position at time t\n";
    std::cout << "  v(t)     = velocity at time t\n\n";

    std::cout << "Energy definition:\n";
    std::cout << "  E(t) = 0.5 * ( x(t)^2 + v(t)^2 )\n";
    std::cout << "For the ideal continuous oscillator, E(t) should be constant.\n";
    std::cout << "With Explicit Euler, small energy drift is expected.\n\n";

// Table header
    std::cout << std::left
              << std::setw(8) << "n"
              << std::setw(10) << "h"
              << std::setw(14) << "t (time)"
              << std::setw(16) << "x(t) (pos)"
              << std::setw(16) << "v(t) (vel)"
              << std::setw(20) << "E(t)=0.5(x^2+v^2)"
              << "\n";

    std::cout << std::string(84, '-') << "\n";

// Helper to print one row (i = index in stored solution arrays)
    auto print_row = [&](std::size_t i) {
        const std::size_t n = i;              // if you store every step, i == n
        const double t = sol.t[i];
        const double x = sol.y[i][0];
        const double v = sol.y[i][1];
        const double E = 0.5 * (x * x + v * v);

        std::cout << std::setw(8) << n
                  << std::setw(10) << h
                  << std::setw(14) << t
                  << std::setw(16) << x
                  << std::setw(16) << v
                  << std::setw(20) << E
                  << "\n";
    };

// Print a readable subset: first, then every K steps, then last.
// If your solution stores EVERY step, this works directly.
// If you only stored some steps, n=i is still a valid index label for the printed rows.
    print_row(0);

    for (std::size_t i = print_every_steps; i + print_every_steps < sol.t.size(); i += print_every_steps) {
        print_row(i);
    }

    if (sol.t.size() > 1) print_row(sol.t.size() - 1);
    return 0;
}