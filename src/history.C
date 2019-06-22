#include "../include/history.h"

std::vector<int64_t> history::indices_for(int64_t index) const {
    std::vector<int64_t> indices(_dims);
    for (int64_t i = 0; i < _dims; ++i) {
        indices[i] = index % _shape[i];
        index = index / _shape[i];
    }

    return indices;
}

void history::add(history const& other, double c1) {
    for (int64_t i = 0; i < _size; ++i)
        histograms[i]->Add(other[i], c1);
}

void history::operator+=(history const& other) { this->add(other, 1); }
void history::operator-=(history const& other) { this->add(other, -1); }

void history::scale(double c1) {
    for (auto const& hist : histograms)
        hist->Scale(c1);
}

void history::operator*=(double c1) { this->scale(c1); }
void history::operator/=(double c1) { this->scale(1. / c1); }

void history::multiply(history const& other) {
    /* assume self, other have equal shapes */
    for (int64_t j = 0; j < _size; ++j) {
        auto count = other[j]->GetBinContent(1);
        histograms[j]->Scale(count);
    }
}

void history::divide(history const& other) {
    /* assume self, other have equal shapes */
    for (int64_t j = 0; j < _size; ++j) {
        auto count = other[j]->GetBinContent(1);
        auto scale = count != 0 ? 1. / count : 0;
        histograms[j]->Scale(scale);
    }
}

void history::operator*=(history const& other) { this->multiply(other); }
void history::operator/=(history const& other) { this->divide(other); }

void history::divide(TH1* const other) {
    for (auto const& hist : histograms)
        hist->Divide(other);
}

TH1F*& history::operator[](int64_t index) {
    return histograms[index];
}

TH1F* const& history::operator[](int64_t index) const {
    return histograms[index];
}

TH1F* history::sum(std::vector<int64_t> const& indices, int64_t axis) const {
    using namespace std::literals::string_literals;

    std::vector<int64_t> output = indices;
    output.erase(std::next(std::begin(output), axis));

    std::string full_tag = _tag + "_sum"s + std::to_string(axis);
    for (auto const& index : output)
        full_tag = full_tag + "_"s + std::to_string(index);

    return sum_impl(full_tag, indices, axis, 0, _shape[axis]);
}

TH1F* history::sum(std::vector<int64_t> const& indices, int64_t axis,
                   int64_t start, int64_t end) const {
    using namespace std::literals::string_literals;

    std::vector<int64_t> output = indices;
    output.erase(std::next(std::begin(output), axis));

    std::string full_tag = _tag + "_sum"s + std::to_string(axis)
        + "f"s + std::to_string(start) + "t"s + std::to_string(end);
    for (auto const& index : output)
        full_tag = full_tag + "_"s + std::to_string(index);

    return sum_impl(full_tag, indices, axis, start, end);
}

void history::apply(std::function<void(TH1*)> f) {
    for (auto& hist : histograms) { f(hist); }
}

TH1F* history::sum_impl(std::string const& name, std::vector<int64_t> indices,
                        int64_t axis, int64_t start, int64_t end) const {
    auto sum = static_cast<TH1F*>((*this)[indices]->Clone(name.data()));
    sum->Reset("MICES");
    for (int64_t i = start; i < end; ++i) {
        indices[axis] = i;
        sum->Add((*this)[indices]);
    }

    return sum;
}

std::unique_ptr<history> history::sum_impl(int64_t axis) const {
    using namespace std::literals::string_literals;

    std::vector<int64_t> output = _shape;
    output.erase(std::next(std::begin(output), axis));

    auto result = std::make_unique<history>(
        _tag + "_sum"s + std::to_string(axis), _ordinate, bins, output);

    for (int64_t i = 0; i < result->size(); ++i) {
        auto indices = result->indices_for(i);
        indices.insert(std::next(std::begin(indices), axis), 0);
        (*result)[i] = this->sum(indices, axis);
    }

    return result;
}

void history::allocate_histograms() {
    using namespace std::literals::string_literals;

    histograms = std::vector<TH1F*>(_size, 0);
    for (int64_t i = 0; i < _size; ++i) {
        std::string index_string;
        for (auto const& index : indices_for(i))
            index_string = index_string + "_"s + std::to_string(index);

        histograms[i] = bins->book<TH1F>(_tag + index_string,
            ";"s + bins->abscissa() + ";"s + _ordinate);
    }
}
