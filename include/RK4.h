//
// Created by Dimitar Vitliyanov on 19.02.26.
//

#ifndef SENIOR_PROJECT_RK4_H
#define SENIOR_PROJECT_RK4_H
#include "Stepper.h"
#include <vector>

class RK4 :Stepper{
public:
    void step(const ODE& ode,
              double t,
              double h,
              const std::vector<double>& y,
              std::vector<double>& y_next) const override;

};


#endif //SENIOR_PROJECT_RK4_H
