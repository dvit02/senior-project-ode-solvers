#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "Euler.h"
#include "RK4.h"
#include "ExponentialDecayODE.h"
#include "Integrator.h"
#include "Solution.h"
int main()
{
    const double k = 1.0;     // decay rate
    const double y0 = 1.0;    // initial condition y(0)
    const double t0 = 0.0;    // start time
    const double tEnd = 2.0;  // end time

    ExponentialDecayODE ode(k); // create exponential decay ODE
    const State init = { y0 };  // wrap scalar into 1-element state vector
    Integrator integrator;      // fixed-step integrator
    Euler euler;                // Euler stepper
    RK4 rk4;                    // RK4 stepper

    // Exact solution
    auto y_exact = [&](double tt) {   // analytical solution function
        return y0 * std::exp(-k * tt); // y(t) = y0 * e^{-kt}
    };

    // Max absolute error over the computed grid
    auto max_abs_error = [&](const Solution& sol) {
        double m = 0.0; // store maximum error
        for (std::size_t i = 0; i < sol.t.size(); ++i)
        {
            m = std::max(m, std::abs(sol.y[i][0] - y_exact(sol.t[i]))); // compare numerical vs exact
        }
        return m; // return maximum error
    };

    std::cout << std::fixed << std::setprecision(6); // set output formatting

    std::cout << "Convergence test on y' = -k y, k=" << k
              << ", y(0)=" << y0
              << ", t in [" << t0 << "," << tEnd << "]\n\n";

    // table header
    std::cout << std::left
              << std::setw(12) << "h"
              << std::setw(22) << "Euler max err"
              << std::setw(22) << "RK4 max err"
              << std::setw(14) << "rate(E)"
              << std::setw(14) << "rate(RK4)"
              << "\n";

    std::cout << std::string(84, '-') << "\n"; // separator line

    double prevE = 0.0, prevR = 0.0; // previous errors (for convergence rate)

    const double h0 = 0.2; // initial step size
    const int Nlevels = 6; // number of refinement levels

    for (int lev = 0; lev < Nlevels; ++lev)
    {
        const double h = h0 / std::pow(2.0, lev); // halve step size each iteration

        const auto solE = integrator.integrate(ode, euler, t0, init, tEnd, h); // Euler solution
        const auto solR = integrator.integrate(ode, rk4,   t0, init, tEnd, h); // RK4 solution

        const double errE = max_abs_error(solE); // Euler max error
        const double errR = max_abs_error(solR); // RK4 max error

        double rateE = 0.0;
        double rateR = 0.0;

        if (lev > 0) // compute convergence rate after first iteration
        {
            rateE = std::log(prevE / errE) / std::log(2.0); // Euler convergence order
            rateR = std::log(prevR / errR) / std::log(2.0); // RK4 convergence order
        }

        // Print h as fixed; errors as scientific with more digits so they don't round to 0
        std::cout << std::setw(12)
                  << std::fixed << std::setprecision(6) << h
                  << std::setw(22)
                  << std::scientific << std::setprecision(10) << errE
                  << std::setw(22)
                  << std::scientific << std::setprecision(10) << errR
                  << std::setw(14)
                  << std::fixed << std::setprecision(6) << (lev == 0 ? 0.0 : rateE)
                  << std::setw(14)
                  << std::fixed << std::setprecision(6) << (lev == 0 ? 0.0 : rateR)
                  << "\n";

        prevE = errE; // store current error for next rate computation
        prevR = errR;
    }

    const double h_detail = 0.1; // step size for detailed comparison table

    std::cout << "\nDetailed table for h = "
              << std::fixed << std::setprecision(6)
              << h_detail << "\n";

    // detailed table header
    std::cout << std::left
              << std::setw(10) << "t"
              << std::setw(16) << "y (Euler)"
              << std::setw(16) << "y (RK4)"
              << std::setw(16) << "y (exact)"
              << std::setw(18) << "|E-exact|"
              << std::setw(18) << "|RK4-exact|"
              << "\n";

    std::cout << std::string(96, '-') << "\n"; // separator

    const auto solE = integrator.integrate(ode, euler, t0, init, tEnd, h_detail); // Euler detailed solution
    const auto solR = integrator.integrate(ode, rk4,   t0, init, tEnd, h_detail); // RK4 detailed solution

    const std::size_t n = std::min(solE.t.size(), solR.t.size()); // ensure same size

    for (std::size_t i = 0; i < n; ++i)
    {
        const double t = solE.t[i]; // current time
        const double ye = y_exact(t); // exact solution

        const double errEi = std::abs(solE.y[i][0] - ye); // Euler pointwise error
        const double errRi = std::abs(solR.y[i][0] - ye); // RK4 pointwise error

        // Keep state values fixed, show errors in scientific so they don't become 0.000000
        std::cout << std::setw(10)
                  << std::fixed << std::setprecision(6) << t
                  << std::setw(16)
                  << std::fixed << std::setprecision(6) << solE.y[i][0]
                  << std::setw(16)
                  << std::fixed << std::setprecision(6) << solR.y[i][0]
                  << std::setw(16)
                  << std::fixed << std::setprecision(6) << ye
                  << std::setw(18)
                  << std::scientific << std::setprecision(10) << errEi
                  << std::setw(18)
                  << std::scientific << std::setprecision(10) << errRi
                  << "\n";
    }

    return 0; // program finished
}
