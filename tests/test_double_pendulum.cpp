//
// Created by Dimitar Vitliyanov on 8.04.26.
//
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "DoublePendulum.h"
#include "RK4.h"
#include "Integrator.h"

namespace {
    bool isFiniteVector(const State& y) {
        for (double v : y) {
            if (!std::isfinite(v)) return false;
        }
        return true;
    }

    double energyDrift(const DoublePendulum& model, const Solution& sol) {
        assert(!sol.y.empty());
        const double E0 = model.energy(sol.y.front());
        const double E1 = model.energy(sol.y.back());
        return std::abs(E1 - E0);
    }
}

int main() {
    // 1. Basic model sanity
    DoublePendulum model(1.0, 1.0, 1.0, 1.0, 9.81);

    assert(model.dim() == 4);

    State y0 = {1.0, 0.5, 0.0, 0.0}; // theta1, theta2, omega1, omega2
    State dydt(model.dim(), 0.0);

    model.rhs(0.0, y0, dydt);

    assert(dydt.size() == 4);
    assert(isFiniteVector(dydt));

    const double E0 = model.energy(y0);
    assert(std::isfinite(E0));

    // -------------------------------------------------------------------------
    // 2. Short-time RK4 refinement sanity:
    //    smaller h should not have worse final energy drift than larger h
    // -------------------------------------------------------------------------
    RK4 rk4;
    Integrator integrator;

    const double t0 = 0.0;
    const double tEnd = 2.0;

    const double h_coarse = 1e-2;
    const double h_fine   = 5e-3;

    Solution solCoarse = integrator.integrate(model, rk4, t0, y0, tEnd, h_coarse);
    Solution solFine   = integrator.integrate(model, rk4, t0, y0, tEnd, h_fine);

    assert(!solCoarse.y.empty());
    assert(!solFine.y.empty());

    const double driftCoarse = energyDrift(model, solCoarse);
    const double driftFine   = energyDrift(model, solFine);

    assert(std::isfinite(driftCoarse));
    assert(std::isfinite(driftFine));

    // Allow equality, but refined RK4 should not be worse here.
    assert(driftFine <= driftCoarse + 1e-12);

    std::cout << "test_double_pendulum passed\n";
    std::cout << "Coarse RK4 drift = " << driftCoarse << "\n";
    std::cout << "Fine   RK4 drift = " << driftFine << "\n";

    return 0;
}