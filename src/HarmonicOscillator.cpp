//
// Created by Dimitar Vitliyanov on 23.02.26.
//

#include "HarmonicOscillator.h"

HarmonicOscillator::HarmonicOscillator(double omega)
        : omega_(omega) {}

std::size_t HarmonicOscillator::dim() const {
    return 2;
}

void HarmonicOscillator::rhs(double /*t*/, const State& y, State& dydt) const {
    if (dydt.size() != 2) dydt.resize(2);

    const double x = y[0];
    const double v = y[1];

    dydt[0] = v;
    dydt[1] = -(omega_ * omega_) * x;
}