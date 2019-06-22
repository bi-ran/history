#ifndef MULTIVAL_H
#define MULTIVAL_H

#include <iterator>
#include <numeric>
#include <type_traits>
#include <vector>

#include "interval.h"

class multival {
  public:
    template <typename... T>
    multival(T const&... intervals)
            : _dims(sizeof...(T)) {
        extract(intervals...);
        _size = std::accumulate(std::begin(_shape), std::end(_shape),
                                1, std::multiplies<int64_t>());
    }

    multival(multival const& other) = default;
    multival& operator=(multival const& other) = default;
    ~multival() = default;

    template <template <typename...> class T>
    std::vector<int64_t> indices_for(T<double> const& values) const {
        std::vector<int64_t> indices;
        auto v = std::begin(values);
        for (auto const& dim : _intervals) {
            indices.push_back(dim.index_for(*v));
            std::advance(v, 1);
        }

        return indices;
    }

    template <template <typename...> class T, typename U>
    typename std::enable_if<std::is_integral<U>::value, int64_t>::type
    index_for(T<U> const& indices) const {
        int64_t index = 0;
        int64_t block = 1;
        auto x = std::begin(indices);
        for (auto const& axis : _shape) {
            index = index + (*x) * block;
            block = block * axis;
            std::advance(x, 1);
        }

        return index;
    }

    template <template <typename...> class T, typename U>
    typename std::enable_if<std::is_floating_point<U>::value, int64_t>::type
    index_for(T<U> const& values) const {
        return index_for(indices_for(values)); }

    std::vector<int64_t> const& shape() const { return _shape; }
    int64_t const& dims() const { return _dims; }
    int64_t const& size() const { return _size; }

  private:
    template <typename... T>
    void extract(T const&... args) {
        (void) (int [sizeof...(T)]) { (_intervals.emplace_back(args), 0)... };
        for (auto const& dim : _intervals) { _shape.push_back(dim.size()); }
    }

    std::vector<int64_t> _shape;
    int64_t _dims;
    int64_t _size;

    std::vector<interval> _intervals;
};

#endif /* MULTIVAL_H */
