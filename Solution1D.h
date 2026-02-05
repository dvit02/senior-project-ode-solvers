//
// Created by Dimitar Vitliyanov on 5.02.26.
//

#ifndef SENIOR_PROJECT_SOLUTION1D_H
#define SENIOR_PROJECT_SOLUTION1D_H
#include <vector>

struct Solution1D {
    std::vector<double> t;   // time points
    std::vector<double> y;   // solution values

    std::size_t size() const {
        return t.size();
    }
};




#endif //SENIOR_PROJECT_SOLUTION1D_H
