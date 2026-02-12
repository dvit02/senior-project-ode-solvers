//
// Created by Dimitar Vitliyanov on 10.02.26.
//

#ifndef SENIOR_PROJECT_STATE_H
#define SENIOR_PROJECT_STATE_H
#include <vector>
#include <cstddef>

// State of an ODE system: y(t) in R^d
using State = std::vector<double>;

// Helper: construct a state of dimension d, optionally with an initial value
inline State make_state(std::size_t d, double value = 0.0) {
    return State(d, value);
}
#endif
