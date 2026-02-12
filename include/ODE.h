//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_ODE_H
#define SENIOR_PROJECT_ODE_H

#include "State.h"
#include <cstddef>
// Interface for an ODE system: y'(t) = f(t, y(t)),  y in R^d
class ODE {
public:
    virtual ~ODE() = default;

    // Dimension d of the state vector
    virtual std::size_t dim() const = 0;

    // Compute rhs: dydt = f(t, y)
    // Requirements:
    //  - y.size() == dim()
    //  - dydt will be resized to dim() (or assumed already correct)
    virtual void rhs(double t, const State& y, State& dydt) const = 0;
};
#endif //SENIOR_PROJECT_ODE_H
