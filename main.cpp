#include <iostream>
#include <iomanip>

#include "Newton1D.h"
#include "TestQuadraticEq.h"

int main() {
    TestQuadraticEq eq;

    Newton1D solver({1e-12, 50});
    auto res = solver.runNewton(eq, 1.0);

    std::cout << std::setprecision(16);
    std::cout << "Converged: " << (res.converged ? "yes" : "no") << "\n";
    std::cout << "Iterations: " << res.iterations << "\n";
    std::cout << "Root: " << res.root << "\n";
    std::cout << "|f(root)|: " << res.finalAbsF << "\n";
    return 0;
}
