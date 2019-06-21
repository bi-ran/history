#include "../include/interval.h"

#include "TH1F.h"
#include "TH1D.h"

#include <iterator>
#include <numeric>

interval::interval(int64_t number, double min, double max,
                   std::string const& abscissa)
        : _size(number),
          _edges(std::vector<double>(number + 1)),
          _abscissa(abscissa) {
    std::iota(std::begin(_edges), std::end(_edges), 0);
    double interval = (max - min) / number;
    for (auto& edge : _edges)
        edge = min + edge * interval;
}

interval::interval(int64_t number, double min, double max)
        : interval(number, min, max, std::string()) { }

int64_t interval::index_for(double value) const {
    int64_t index = _size;
    for (auto edge : _edges)
        if (value < edge)
            --index;

    return index;
}

template <typename T>
T* interval::book(std::string const& name, std::string const& title) {
    return new T(name.data(), title.data(), _size, _edges.data());
}

/* explicit instantiations */
template TH1F* interval::book<TH1F>(std::string const&, std::string const&);
template TH1D* interval::book<TH1D>(std::string const&, std::string const&);
