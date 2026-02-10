//
// Created by Dimitar Vitliyanov on 29.01.26.
//

#include "../include/ScalarFunction.h"
#include <cmath>
using namespace std;
double ScalarFunction::derivative(double x) const {
    const double epsilon = 1e-6 * (1.0 + abs(x));// this is the h = 10^-6(1+|x|)
    return (value(x + epsilon) - value(x - epsilon)) / (2.0 * epsilon);//deffinitio of a derivative
}
