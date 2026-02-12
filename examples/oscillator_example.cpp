//
// Created by Dimitar Vitliyanov on 10.02.26.
//
#include <iostream>
#include <iomanip>
#include <cmath>

#include "ODE.h"
#include "Integrator.h"
#include "Euler.h"

// Harmonic oscillator: x' = v, v' = -x
class HarmonicOscillator : public ODE {
public:
    std::size_t dim() const override { return 2; }

    void rhs(double /*t*/, const State& y, State& dydt) const override {
        if (dydt.size() != 2) dydt.resize(2);
        // y[0] = x, y[1] = v
        dydt[0] = y[1];
        dydt[1] = -y[0];
    }
};

int main() {
    HarmonicOscillator ode;
    Euler stepper;
    Integrator integrator;

    State y0 = {1.0, 0.0}; // x(0)=1, v(0)=0
    Solution sol = integrator.integrate(ode, stepper, 0.0, y0, 10.0, 0.01);

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\nHarmonic Oscillator (Euler)\n";
    std::cout << "System: x' = v,  v' = -x\n";
    std::cout << "Initial: x(0)=1, v(0)=0\n\n";

    std::cout << std::setw(10) << "t"
              << std::setw(14) << "x"
              << std::setw(14) << "v"
              << std::setw(18) << "E=0.5(x^2+v^2)"
              << "\n";

    std::cout << std::string(10 + 14 + 14 + 18, '-') << "\n";

    for (std::size_t i = 0; i < sol.t.size(); i += 200) {
        const double t = sol.t[i];
        const double x = sol.y[i][0];
        const double v = sol.y[i][1];
        const double E = 0.5 * (x * x + v * v);

        std::cout << std::setw(10) << t
                  << std::setw(14) << x
                  << std::setw(14) << v
                  << std::setw(18) << E
                  << "\n";
    }

    std::cout << "\n";
    return 0;
}