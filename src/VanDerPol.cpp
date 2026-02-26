//
// Created by Dimitar Vitliyanov on 25.02.26.
//
//
// Created by Dimitar Vitliyanov on 25.02.26.
//

#include "VanDerPol.h"

VanDerPol::VanDerPol(double mu)
        : mu_(mu)
{}

std::size_t VanDerPol::dim() const
{
    return 2;
}

void VanDerPol::rhs(double /*t*/, const State& y, State& dydt) const
{
    const double x = y[0];
    const double v = y[1];

    dydt.resize(2);
    dydt[0] = v;
    dydt[1] = mu_ * (1.0 - x * x) * v - x;
}