//
// Created by Dimitar Vitliyanov on 29.01.26.
//

#include "../../examples/TestQuadraticEq.h"

double TestQuadraticEq::value(double x) const {
    return x * x - 2.0;
}

double TestQuadraticEq::derivative(double x) const {
    return 2.0 * x;
}
