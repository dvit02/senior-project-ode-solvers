//
// Created by Dimitar Vitliyanov on 25.02.26.
//

#ifndef SENIOR_PROJECT_VANDERPOL_H
#define SENIOR_PROJECT_VANDERPOL_H
#include "ODE.h"

class VanDerPol : public ODE {
public:
    explicit VanDerPol(double mu);

    std::size_t dim() const override;
    void rhs(double t, const State& y, State& dydt) const override;

private:
    double mu_;
};

#endif //SENIOR_PROJECT_VANDERPOL_H
