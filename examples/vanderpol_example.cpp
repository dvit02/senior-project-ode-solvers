//
// Created by Dimitar Vitliyanov on 25.02.26.
//
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Integrator.h"
#include "RK4.h"
#include "VanDerPol.h"

// Export fixed-step solution
static void export_fixed(const char* filename, const Solution& sol)
{
    std::ofstream out(filename);
    out << "t,x,v\n";

    for (std::size_t i = 0; i < sol.t.size(); ++i)
    {
        out << sol.t[i] << ","
            << sol.y[i][0] << ","
            << sol.y[i][1] << "\n";
    }
}

// Export adaptive solution (with step size info)
static void export_adaptive(const char* filename,
                            const AdaptiveResult& result)
{
    std::ofstream out(filename);
    out << "t,x,v\n";

    const Solution& sol = result.solution;

    for (std::size_t i = 0; i < sol.t.size(); ++i)
    {
        out << sol.t[i] << ","
            << sol.y[i][0] << ","
            << sol.y[i][1] << "\n";
    }
}

int main()
{
    const double mu = 10.0;

    VanDerPol ode(mu);
    RK4 stepper;
    Integrator integrator;

    const double t0   = 0.0;
    const double tEnd = 20.0;

    State y0 = {2.0, 0.0};

    // Fixed-step integration
    const double h_fixed = 0.01;

    Solution fixed =
            integrator.integrate(ode, stepper, t0, y0, tEnd, h_fixed);

    export_fixed("vdp_fixed.csv", fixed);

    // Adaptive integration
    const double h0   = 0.05;
    const double rtol = 1e-6;
    const double atol = 1e-9;
    const double h_min = 1e-6;
    const double h_max = 0.1;

    AdaptiveResult adaptive =
            integrator.integrateAdaptiveRK4(
                    ode, stepper,
                    t0, y0, tEnd,
                    h0, rtol, atol,
                    h_min, h_max
            );

    export_adaptive("vdp_adaptive.csv", adaptive);

    // Diagnostics
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Van der Pol oscillator (mu = " << mu << ")\n\n";

    std::cout << "Fixed step size: h = " << h_fixed << "\n";
    std::cout << "Fixed total steps: "
              << fixed.t.size() - 1 << "\n\n";

    std::cout << "Adaptive integration:\n";
    std::cout << "Accepted steps: "
              << adaptive.accepted_steps << "\n";
    std::cout << "Rejected steps: "
              << adaptive.rejected_steps << "\n";
    std::cout << "Minimum h used: "
              << adaptive.h_min_used << "\n";
    std::cout << "Maximum h used: "
              << adaptive.h_max_used << "\n";

    return 0;
}