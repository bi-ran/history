#include "../include/multival.h"

#include "TH2.h"
#include "TH3.h"

using namespace std::literals::string_literals;

std::vector<int64_t> multival::indices_for(int64_t index) const {
    std::vector<int64_t> indices(_dims);
    for (int64_t i = 0; i < _dims; ++i) {
        indices[i] = index % _shape[i];
        index = index / _shape[i];
    }

    return indices;
}

/* template specialisations */
template <>
TH2F* multival::book<TH2F>(int64_t, std::string const& name,
                           std::string const& ordinate) const {
    auto title = ";"s + _intervals[0].abscissa() + ";"s
        + _intervals[1].abscissa() + ";"s + ordinate;

    return new TH2F(name.data(), title.data(),
        _intervals[0].size(), _intervals[0].edges(),
        _intervals[1].size(), _intervals[1].edges());
}

template <>
TH3F* multival::book<TH3F>(int64_t, std::string const& name,
                           std::string const&) const {
    auto title = ";"s + _intervals[0].abscissa() + ";"s
        + _intervals[1].abscissa() + ";"s + _intervals[2].abscissa();

    return new TH3F(name.data(), title.data(),
        _intervals[0].size(), _intervals[0].edges(),
        _intervals[1].size(), _intervals[1].edges(),
        _intervals[2].size(), _intervals[2].edges());
}

/* explicit instantiations */
template TH2F*
multival::book<TH2F>(int64_t, std::string const&, std::string const&) const;
template TH3F*
multival::book<TH3F>(int64_t, std::string const&, std::string const&) const;
