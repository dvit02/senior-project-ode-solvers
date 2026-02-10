//
// Created by Dimitar Vitliyanov on 29.01.26.
//

#ifndef SENIOR_PROJECT_NEWTON1D_H
#define SENIOR_PROJECT_NEWTON1D_H

#include "ScalarFunction.h"

struct Newton1DOptions { // we need data container, hence we use struct
    double tol = 1e-10;  // stopping tolerance on |f(x)|; how close to zero is from f(x)
    int maxIter = 50;    // maximum number of iterations
};

struct Newton1DResult {
    double root = 0.0; //final x reached by newton
    int iterations = 0;
    bool converged = false;
    double finalAbsF = 0.0;
    // what it returns after running
};

class Newton1D {
public:
    explicit Newton1D(Newton1DOptions options = Newton1DOptions()); //explicit forces construction to be intentional, not accidental.


    // Solve f(x)=0 starting from x0
    Newton1DResult runNewton (const ScalarFunction& f, double x0) const;

private:
    Newton1DOptions opt_;
//It stores the options inside the solver, so you don’t need to pass tolerance/maxIter every time you call runNewton .
};


#endif //SENIOR_PROJECT_NEWTON1D_H
