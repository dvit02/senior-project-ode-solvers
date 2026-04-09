//
// Created by Dimitar Vitliyanov on 8.04.26.
//

#include "DoublePendulum.h"
#include <cmath>

DoublePendulum::DoublePendulum(double m1, double m2, double L1, double L2, double g)
        : m1_(m1), m2_(m2), L1_(L1), L2_(L2), g_(g)
{}

std::size_t DoublePendulum::dim() const
{
    return 4; // theta1, theta2, omega1, omega2
}

void DoublePendulum::rhs(double /*t*/, const State& y, State& dydt) const
{
    const double th1 = y[0]; // angle of first bob from vertical
    const double th2 = y[1]; // angle of second bob from vertical
    const double w1  = y[2]; // angular velocity of first bob
    const double w2  = y[3]; // angular velocity of second bob

    const double d   = th1 - th2;         // angle difference shows up a lot
    const double s   = std::sin(d);       // precompute to avoid recomputing
    const double c   = std::cos(d);

    const double M   = m1_ + m2_;         // total mass shortcut
    const double denom_common = 2.0 * m1_ + m2_ - m2_ * std::cos(2.0 * th1 - 2.0 * th2);
    // shared denominator from the Lagrangian-derived form

    // standard closed-form equations of motion (see Diego Assencio derivation)
    const double num1 =
            -g_ * (2.0 * m1_ + m2_) * std::sin(th1)
            - m2_ * g_ * std::sin(th1 - 2.0 * th2)
            - 2.0 * s * m2_ * (w2 * w2 * L2_ + w1 * w1 * L1_ * c);

    const double num2 =
            2.0 * s * (w1 * w1 * L1_ * M
                       + g_ * M * std::cos(th1)
                       + w2 * w2 * L2_ * m2_ * c);

    dydt.resize(4);
    dydt[0] = w1;                          // theta1 dot
    dydt[1] = w2;                          // theta2 dot
    dydt[2] = num1 / (L1_ * denom_common); // omega1 dot
    dydt[3] = num2 / (L2_ * denom_common); // omega2 dot
}

double DoublePendulum::energy(const State& y) const
{
    const double th1 = y[0];
    const double th2 = y[1];
    const double w1  = y[2];
    const double w2  = y[3];

    // kinetic energy of both bobs, second bob uses coupled velocity
    const double T =
            0.5 * m1_ * L1_ * L1_ * w1 * w1
            + 0.5 * m2_ * (L1_ * L1_ * w1 * w1
                           + L2_ * L2_ * w2 * w2
                           + 2.0 * L1_ * L2_ * w1 * w2 * std::cos(th1 - th2));

    // potential energy measured from pivot, so both terms are negative when hanging
    const double V =
            -(m1_ + m2_) * g_ * L1_ * std::cos(th1)
            - m2_ * g_ * L2_ * std::cos(th2);

    return T + V;
}