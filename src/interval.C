#include "../include/interval.h"

#include "TH1F.h"
#include "TH2F.h"

#include <iterator>
#include <numeric>

interval::interval(std::string const& abscissa, int64_t number,
                   double min, double max)
        : _abscissa(abscissa),
          _size(number),
          _edges(std::vector<double>(number + 1)) {
    std::iota(std::begin(_edges), std::end(_edges), 0);
    double interval = (max - min) / number;
    for (auto& edge : _edges)
        edge = min + edge * interval;
}

interval::interval(int64_t number, double min, double max)
        : interval(std::string(), number, min, max) { }

int64_t interval::index_for(double value) const {
    int64_t index = _size;
    for (auto edge : _edges)
        if (value < edge)
            --index;

    return index;
}

std::array<double, 2> interval::edges(int64_t index) const {
    return { _edges[index], _edges[index + 1] };
}

/* template specialisations */

template <>
TH1F* interval::book<TH1F>(std::string const& name,
                           std::string const& title) const {
    return new TH1F(name.data(), title.data(), _size, _edges.data());
}

template <>
TH2F* interval::book<TH2F>(std::string const& name,
                           std::string const& title) const {
    return new TH2F(name.data(), title.data(), _size, _edges.data(),
                    _size, _edges.data());
}

/* explicit instantiations */
template TH1F*
interval::book<TH1F>(std::string const&, std::string const&) const;
template TH2F*
interval::book<TH2F>(std::string const&, std::string const&) const;
