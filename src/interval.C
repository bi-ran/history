#include "../include/interval.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <iterator>
#include <numeric>

using namespace std::literals::string_literals;

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

/* template specialisations */
template <>
TH1F* interval::book<TH1F>(int64_t, std::string const& name,
                           std::string const& ordinate) const {
    auto title = ";"s + _abscissa + ";"s + ordinate;
    return new TH1F(name.data(), title.data(), _size, _edges.data());
}

template <>
TH2F* interval::book<TH2F>(int64_t, std::string const& name,
                           std::string const& ordinate) const {
    auto title = ";"s + _abscissa + ";"s + _abscissa + ";"s + ordinate;
    return new TH2F(name.data(), title.data(), _size, _edges.data(),
        _size, _edges.data());
}

template <>
TH3F* interval::book<TH3F>(int64_t, std::string const& name,
                           std::string const&) const {
    auto title = ";"s + _abscissa + ";"s + _abscissa + ";"s + _abscissa;
    return new TH3F(name.data(), title.data(), _size, _edges.data(),
        _size, _edges.data(), _size, _edges.data());
}

template <>
TH1F* interval::book<TH1F, 2>(int64_t, std::string const& name,
                              std::string const& ordinate,
                              std::array<int64_t, 2> const& offset) const {
    auto title = ";"s + _abscissa + ";"s + ordinate;
    return new TH1F(name.data(), title.data(), _size - offset[0] - offset[1],
                    _edges.data() + offset[0]);
}

template <>
TH2F* interval::book<TH2F, 4>(int64_t, std::string const& name,
                              std::string const& ordinate,
                              std::array<int64_t, 4> const& offset) const {
    auto title = ";"s + _abscissa + ";"s + _abscissa + ";"s + ordinate;
    return new TH2F(name.data(), title.data(), _size - offset[0] - offset[1],
                    _edges.data() + offset[0], _size - offset[2] - offset[3],
                    _edges.data() + offset[2]);
}

template <>
TH3F* interval::book<TH3F, 6>(int64_t, std::string const& name,
                              std::string const&,
                              std::array<int64_t, 6> const& offset) const {
    auto title = ";"s + _abscissa + ";"s + _abscissa + ";"s + _abscissa;
    return new TH3F(name.data(), title.data(), _size - offset[0] - offset[1],
                    _edges.data() + offset[0], _size - offset[2] - offset[3],
                    _edges.data() + offset[2], _size - offset[4] - offset[5],
                    _edges.data() + offset[4]);
}

/* explicit instantiations */
template TH1F*
interval::book<TH1F>(int64_t, std::string const&, std::string const&) const;
template TH2F*
interval::book<TH2F>(int64_t, std::string const&, std::string const&) const;
template TH3F*
interval::book<TH3F>(int64_t, std::string const&, std::string const&) const;
template TH1F*
interval::book<TH1F, 2>(int64_t, std::string const&, std::string const&,
                        std::array<int64_t, 2> const&) const;
template TH2F*
interval::book<TH2F, 4>(int64_t, std::string const&, std::string const&,
                        std::array<int64_t, 4> const&) const;
template TH3F*
interval::book<TH3F, 6>(int64_t, std::string const&, std::string const&,
                        std::array<int64_t, 6> const&) const;
