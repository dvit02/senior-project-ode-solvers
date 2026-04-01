//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_EXPONENTIALDECAYODE_H
#define SENIOR_PROJECT_EXPONENTIALDECAYODE_H
#include "ODE.h"

class ExponentialDecayODE: public ODE {
public:
    explicit ExponentialDecayODE(double k);

    std::size_t dim() const override;

    void rhs(double t, const State& y, State& dydt) const override;
private:
    double k_;
};


#endif //SENIOR_PROJECT_EXPONENTIALDECAYODE_H
