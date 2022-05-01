#include "tools.hpp"
#include <algorithm>

//#include <iostream>

namespace utils {
    float clip(float n, float lower, float upper) {
        return std::max(lower, std::min(n, upper));
    }

    int IX(const int x, const int y, const int size) {
         return clip(x, 0, size - 1) + (clip(y, 0, size - 1) * size);
    }
}