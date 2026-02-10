//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_ODE1D_H
#define SENIOR_PROJECT_ODE1D_H


class ODE1D {
public:
    virtual ~ODE1D() = default;

    // returns f(t, y) = dy/dt
    virtual double rhs(double t, double y) const = 0;

};


#endif //SENIOR_PROJECT_ODE1D_H
