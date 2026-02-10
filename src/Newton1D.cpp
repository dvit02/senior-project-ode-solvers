//
// Created by Dimitar Vitliyanov on 29.01.26.
//

#include "../include/Newton1D.h"
#include <cmath>
#include <stdexcept>

Newton1D::Newton1D(Newton1DOptions options) : opt_(options) {}

Newton1DResult Newton1D::runNewton(const ScalarFunction &f, double x0) const {
    double x = x0;

    for (int it = 1; it <= opt_.maxIter; ++it) {
        const double fx  = f.value(x);
        const double dfx = f.derivative(x);

        if (!std::isfinite(fx) || !std::isfinite(dfx)) {
            throw std::runtime_error("Newton1D: non-finite f(x) or f'(x).");
        }
        if (std::abs(dfx) < 1e-14) {
            throw std::runtime_error("Newton1D: derivative too small (division by ~0).");
        }

        // Newton update: x_{n+1} = x_n - f(x_n)/f'(x_n)
        const double x_new = x - fx / dfx;

        // Convergence check (primary): |f(x)| small
        if (std::abs(fx) < opt_.tol) {
            return Newton1DResult{x, it, true, std::abs(fx)};
        }

        x = x_new;
    }

    return Newton1DResult{x, opt_.maxIter, false, std::abs(f.value(x))};
}
