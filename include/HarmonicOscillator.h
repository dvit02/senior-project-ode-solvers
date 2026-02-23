//
// Created by Dimitar Vitliyanov on 23.02.26.
//

#ifndef SENIOR_PROJECT_HARMONICOSCILLATOR_H
#define SENIOR_PROJECT_HARMONICOSCILLATOR_H
#include "ODE.h"

class HarmonicOscillator: public ODE{
    // Harmonic oscillator:
// x' = v
// v' = -omega^2 * x
public:
    explicit HarmonicOscillator(double omega = 1.0);

    std::size_t dim() const override;
    void rhs(double t, const State& y, State& dydt) const override;

private:
    double omega_;

};


#endif //SENIOR_PROJECT_HARMONICOSCILLATOR_H
