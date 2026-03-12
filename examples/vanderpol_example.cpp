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
    std::ofstream out(filename); // open output file
    out << "t,x,v\n"; // CSV header

    for (std::size_t i = 0; i < sol.t.size(); ++i) // loop through all stored solution points
    {
        out << sol.t[i] << "," // write time value
            << sol.y[i][0] << "," // write x component of state
            << sol.y[i][1] << "\n"; // write v component of state
    }
}

// Export adaptive solution (with step size info)
static void export_adaptive(const char* filename,
                            const AdaptiveResult& result)
{
    std::ofstream out(filename); // open output file
    out << "t,x,v\n"; // CSV header

    const Solution& sol = result.solution; // reference to solution stored inside AdaptiveResult

    for (std::size_t i = 0; i < sol.t.size(); ++i) // loop through adaptive solution
    {
        out << sol.t[i] << "," // write time
            << sol.y[i][0] << "," // write x
            << sol.y[i][1] << "\n"; // write v
    }
}

int main()
{
    const double mu = 10.0; // Van der Pol parameter

    VanDerPol ode(mu); // construct Van der Pol ODE system
    RK4 stepper; // Runge-Kutta 4 stepper
    Integrator integrator; // integrator object

    //time domain
    const double t0   = 0.0;
    const double tEnd = 20.0;

    State y0 = {2.0, 0.0}; //initial condition

    const double h_fixed = 0.01;  // Fixed total steps = 20/0.01 = 2000

    Solution fixed =
            integrator.integrate(ode, stepper, t0, y0, tEnd, h_fixed); // run fixed-step integration

    export_fixed("vdp_fixed.csv", fixed); // export fixed-step results

    // Adaptive integration
    const double h0   = 0.05; // initial step size - adaptive start
    const double rtol = 1e-6; // relative tolerance
    const double atol = 1e-9; // absolute tolerance
    //step size limits
    const double h_min = 1e-6;
    const double h_max = 0.1;

    AdaptiveResult adaptive =
            integrator.integrateAdaptiveRK4(
                    ode, stepper,
                    t0, y0, tEnd,
                    h0, rtol, atol,
                    h_min, h_max
            ); // run adaptive RK4 integrator

    export_adaptive("vdp_adaptive.csv", adaptive); // export adaptive results

    // Diagnostics
    std::cout << std::fixed << std::setprecision(6); // format output

    std::cout << "Van der Pol oscillator (mu = " << mu << ")\n\n";

    std::cout << "Fixed step size: h = " << h_fixed << "\n";
    std::cout << "Fixed total steps: "
              << fixed.t.size() - 1 << "\n\n"; // number of steps performed in fixed integration

    std::cout << "Adaptive integration:\n";
    std::cout << "Accepted steps: "
              << adaptive.accepted_steps << "\n"; // successful steps

    std::cout << "Rejected steps: "
              << adaptive.rejected_steps << "\n"; // steps rejected due to large error

    std::cout << "Minimum h used: "
              << adaptive.h_min_used << "\n"; // smallest step size used during adaptive integration

    std::cout << "Maximum h used: "
              << adaptive.h_max_used << "\n"; // largest step size used during adaptive integration

    return 0;}