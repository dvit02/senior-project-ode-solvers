//
// Created by Dimitar Vitliyanov on 5.02.26.
//
#include "ExponentialDecayODE.h"

ExponentialDecayODE::ExponentialDecayODE(double k) : k_(k) {}

std::size_t ExponentialDecayODE::dim() const {
    return 1;   // scalar ODE — one state component
}

void ExponentialDecayODE::rhs(double /*t*/, const State& y, State& dydt) const {
    if (dydt.size() != 1) dydt.resize(1);
    dydt[0] = -k_ * y[0];   // y' = -k * y
}