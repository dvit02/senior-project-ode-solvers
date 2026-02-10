//
// Created by Dimitar Vitliyanov on 29.01.26.
//

#ifndef SENIOR_PROJECT_TESTQUADRATICEQ_H
#define SENIOR_PROJECT_TESTQUADRATICEQ_H
#include "../Senior_Project/include/ScalarFunction.h"
class TestQuadraticEq: public ScalarFunction{
        double value(double x) const override;       // f(x)
        double derivative(double x) const override;  // f'(x) = 2x (analytic)
        bool hasAnalyticDerivative() const override { return true; }



};


#endif //SENIOR_PROJECT_TESTQUADRATICEQ_H
