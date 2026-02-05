//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_EXPONENTIALDECAYODE_H
#define SENIOR_PROJECT_EXPONENTIALDECAYODE_H
#include "ODE1D.h"

class ExponentialDecayODE: public ODE1D {
    public:
        explicit ExponentialDecayODE(double k);

        double rhs(double t, double y) const override;

    private:
        double k_;

};


#endif //SENIOR_PROJECT_EXPONENTIALDECAYODE_H
