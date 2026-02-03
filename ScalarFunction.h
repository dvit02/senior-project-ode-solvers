//
// Created by Dimitar Vitliyanov on 29.01.26.
//

#ifndef SENIOR_PROJECT_SCALARFUNCTION_H
#define SENIOR_PROJECT_SCALARFUNCTION_H


class ScalarFunction {
    public:
        virtual ~ScalarFunction() = default;

        virtual double value(double x) const = 0; //f(x)

        virtual double derivative(double x) const;    // f'(x) - for analytical derivative

        virtual bool hasAnalyticDerivative() const { return false; } // checks if subcalss give analytical derivative
    };

#endif //SENIOR_PROJECT_SCALARFUNCTION_H
