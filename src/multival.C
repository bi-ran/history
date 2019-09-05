#include "../include/multival.h"

#include "TH1.h"
#include "TH2.h"

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
TH1F* multival::book<TH1F>(std::string const& name,
                           std::string const& ordinate) const {
    auto title = ";"s + _intervals[0].abscissa() + ";"s + ordinate;

    return new TH1F(name.data(), title.data(),
        _intervals[0].size(), _intervals[0].edges());
}

template <>
TH2F* multival::book<TH2F>(std::string const& name,
                           std::string const&) const {
    auto title = ";"s + _intervals[0].abscissa()
        + ";"s + _intervals[1].abscissa();

    return new TH2F(name.data(), title.data(),
        _intervals[0].size(), _intervals[0].edges(),
        _intervals[1].size(), _intervals[1].edges());
}

/* explicit instantiations */
template TH1F*
multival::book<TH1F>(std::string const&, std::string const&) const;
template TH2F*
multival::book<TH2F>(std::string const&, std::string const&) const;
