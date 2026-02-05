#include <iostream>
#include <iomanip>
#include <cmath>
#include "Euler1D.h"
#include "ExponentialDecayODE.h"
#include "Integrator1D.h"
#include "Solution1D.h"

int main() {
    // Problem: y' = -k y, y(0) = y0
    const double k  = 1.0;
    const double y0 = 1.0;

    // Integration setup
    const double t0   = 0.0;
    const double tEnd = 2.0;
    const double h    = 0.1;

    ExponentialDecayODE ode(k);
    Euler1D euler;

    double t = t0;
    double y = y0;

    // Table header
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Euler method for y' = -k y,  k=" << k
              << ",  y(0)=" << y0 << ",  h=" << h << "\n\n";

    std::cout << std::left
              << std::setw(10) << "t"
              << std::setw(16) << "y (Euler)"
              << std::setw(16) << "y (exact)"
              << std::setw(16) << "abs error"
              << "\n";

    std::cout << std::string(58, '-') << "\n";

    auto y_exact = [&](double tt) {
        return y0 * std::exp(-k * tt);
    };

    // Print initial row
    {
        const double ye = y_exact(t);
        std::cout << std::setw(10) << t
                  << std::setw(16) << y
                  << std::setw(16) << ye
                  << std::setw(16) << std::abs(y - ye)
                  << "\n";
    }

    auto sol = integrate(euler, ode, t0, y0, tEnd, h);
    for (std::size_t i = 0; i < sol.size(); ++i) {
        const double ye = y_exact(sol.t[i]);

        std::cout << std::setw(10) << sol.t[i]
                  << std::setw(16) << sol.y[i]
                  << std::setw(16) << ye
                  << std::setw(16) << std::abs(sol.y[i] - ye)
                  << "\n";
    }



    return 0;
}
