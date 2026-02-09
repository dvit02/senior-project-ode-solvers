#include <iostream>
#include <iomanip>
#include <cmath>
#include "Euler1D.h"
#include "RK4_1D.h"
#include "ExponentialDecayODE.h"
#include "Integrator1D.h"
#include "Solution1D.h"

int main() {
    // Problem: y' = -k y, y(0) = y0
    const double k  = 1.0;
    const double y0 = 1.0;

    // Integration interval
    const double t0   = 0.0;
    const double tEnd = 2.0;

    ExponentialDecayODE ode(k);
    Euler1D euler;
    RK4_1D  rk4;

    // Exact solution
    auto y_exact = [&](double tt) {
        return y0 * std::exp(-k * tt);
    };

    // Max absolute error over the computed grid
    auto max_abs_error = [&](const Solution1D& sol) {
        double m = 0.0;
        for (std::size_t i = 0; i < sol.size(); ++i) {
            m = std::max(m, std::abs(sol.y[i] - y_exact(sol.t[i])));
        }
        return m;
    };

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Convergence test on y' = -k y,  k=" << k
              << ",  y(0)=" << y0 << ",  t in [" << t0 << "," << tEnd << "]\n\n";

    std::cout << std::left
              << std::setw(12) << "h"
              << std::setw(18) << "Euler max err"
              << std::setw(18) << "RK4 max err"
              << std::setw(14) << "rate(E)"
              << std::setw(14) << "rate(RK4)"
              << "\n";
    std::cout << std::string(76, '-') << "\n";

    double prevE = 0.0, prevR = 0.0;

    const double h0 = 0.2;
    const int Nlevels = 6;

    for (int lev = 0; lev < Nlevels; ++lev) {
        const double h = h0 / std::pow(2.0, lev);

        const auto solE = integrate(euler, ode, t0, y0, tEnd, h);
        const auto solR = integrate(rk4,  ode, t0, y0, tEnd, h);

        const double errE = max_abs_error(solE);
        const double errR = max_abs_error(solR);

        double rateE = 0.0;
        double rateR = 0.0;
        if (lev > 0) {
            rateE = std::log(prevE / errE) / std::log(2.0);
            rateR = std::log(prevR / errR) / std::log(2.0);
        }

        std::cout << std::setw(12) << h
                  << std::setw(18) << errE
                  << std::setw(18) << errR
                  << std::setw(14) << (lev == 0 ? 0.0 : rateE)
                  << std::setw(14) << (lev == 0 ? 0.0 : rateR)
                  << "\n";

        prevE = errE;
        prevR = errR;
    }


    const double h_detail = 0.1;

    std::cout << "\nDetailed table for h = " << h_detail << "\n";

    std::cout << std::left
              << std::setw(10) << "t"
              << std::setw(16) << "y (Euler)"
              << std::setw(16) << "y (RK4)"
              << std::setw(16) << "y (exact)"
              << std::setw(16) << "|E-exact|"
              << std::setw(16) << "|RK4-exact|"
              << "\n";
    std::cout << std::string(90, '-') << "\n";

    const auto solE = integrate(euler, ode, t0, y0, tEnd, h_detail);
    const auto solR = integrate(rk4,  ode, t0, y0, tEnd, h_detail);

    // Same grid because same h (fixed step)
    const std::size_t n = std::min(solE.size(), solR.size());
    for (std::size_t i = 0; i < n; ++i) {
        const double t  = solE.t[i];
        const double ye = y_exact(t);

        std::cout << std::setw(10) << t
                  << std::setw(16) << solE.y[i]
                  << std::setw(16) << solR.y[i]
                  << std::setw(16) << ye
                  << std::setw(16) << std::abs(solE.y[i] - ye)
                  << std::setw(16) << std::abs(solR.y[i] - ye)
                  << "\n";
    }

    return 0;
}
