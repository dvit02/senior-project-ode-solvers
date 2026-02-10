//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#include "../include/ExponentialDecayODE.h"
ExponentialDecayODE::ExponentialDecayODE(double k) : k_(k) {}

double ExponentialDecayODE::rhs(double /*t*/, double y) const {
    return -k_ * y;}
