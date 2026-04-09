//
// Created by Dimitar Vitliyanov on 8.04.26.
//

#ifndef SENIOR_PROJECT_DOUBLEPENDULUM_H
#define SENIOR_PROJECT_DOUBLEPENDULUM_H
#include "ODE.h"

class DoublePendulum : public ODE {
public:
    DoublePendulum(double m1, double m2, double L1, double L2, double g = 9.81);

    std::size_t dim() const override; // 4D state: theta1, theta2, omega1, omega2
    void rhs(double t, const State& y, State& dydt) const override;

    // total mechanical energy (for drift diagnostics)
    double energy(const State& y) const;

    // accessors for the example driver
    double m1() const { return m1_; }
    double m2() const { return m2_; }
    double L1() const { return L1_; }
    double L2() const { return L2_; }
    double g()  const { return g_;  }

private:
    double m1_, m2_; // masses of the two bobs
    double L1_, L2_; // rod lengths
    double g_;       // gravitational acceleration
};

#endif //SENIOR_PROJECT_DOUBLEPENDULUM_H
