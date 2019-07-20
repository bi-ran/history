#include "../include/multival.h"

std::vector<int64_t> multival::indices_for(int64_t index) const {
    std::vector<int64_t> indices(_dims);
    for (int64_t i = 0; i < _dims; ++i) {
        indices[i] = index % _shape[i];
        index = index / _shape[i];
    }

    return indices;
}
