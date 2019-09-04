#include "../include/history.h"

#include "TNamed.h"

using namespace std::literals::string_literals;

history::history(TFile* f, std::string const& tag)
        : _tag(tag) {
    std::string desc = ((TNamed*)f->Get(tag.data()))->GetTitle();
    while (!desc.empty()) {
        desc.erase(0, 1);

        auto pos = desc.find("_"s);
        auto token = desc.substr(0, pos);
        _shape.push_back(std::atoi(token.data()));

        desc.erase(0, pos);
    }

    _dims = _shape.size();
    _size = std::accumulate(std::begin(_shape), std::end(_shape), 1,
                            std::multiplies<int64_t>());

    histograms = std::vector<TH1F*>(_size, nullptr);
    for (int64_t i = 0; i < _size; ++i) {
        auto name = _tag + stub(indices_for(i));
        histograms[i] = (TH1F*)f->Get(name.data());
        histograms[i]->SetName(name.data());
    }
}

history::history(history const& other, std::string const& prefix) 
        : _tag(prefix + "_"s + other._tag),
          _ordinate(other._ordinate),
          _dims(other._dims),
          _size(other._size),
          _shape(other._shape),
          bins(other.bins) {
    for (auto const& hist : other.histograms) {
        auto name = prefix + "_"s + hist->GetName();
        histograms.emplace_back((TH1F*)hist->Clone(name.data()));
        histograms.back()->SetName(name.data());
    }
}

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

void history::operator+=(history const& other) { add(other, 1); }
void history::operator-=(history const& other) { add(other, -1); }

void history::scale(double c1) {
    for (auto const& hist : histograms)
        hist->Scale(c1);
}

void history::operator*=(double c1) { scale(c1); }
void history::operator/=(double c1) { scale(1. / c1); }

void history::multiply(history const& other) {
    /* warning: fails silently! */
    if (!compatible(other)) { return; }

    std::vector<int64_t> axes(_dims - other._dims);
    std::iota(std::begin(axes), std::end(axes), other._dims);

    return multiply(other, axes);
}

void history::divide(history const& other) {
    /* warning: fails silently! */
    if (!compatible(other)) { return; }

    std::vector<int64_t> axes(_dims - other._dims);
    std::iota(std::begin(axes), std::end(axes), other._dims);

    return divide(other, axes);
}

void history::operator*=(history const& other) { multiply(other); }
void history::operator/=(history const& other) { divide(other); }

void history::multiply(TH1* const other) {
    apply([&](TH1* h) { h->Multiply(other); });
}

void history::divide(TH1* const other) {
    apply([&](TH1* h) { h->Divide(other); });
}

void history::operator*=(TH1* const other) { multiply(other); }
void history::operator/=(TH1* const other) { divide(other); }

TH1F*& history::operator[](int64_t index) {
    return histograms[index];
}

TH1F* const& history::operator[](int64_t index) const {
    return histograms[index];
}

TH1F* history::sum(std::vector<int64_t> indices, int64_t axis) const {
    std::vector<int64_t> output = indices;
    output.erase(std::next(std::begin(output), axis));

    auto name = _tag + "_sum"s + std::to_string(axis) + stub(output);
    auto sum = static_cast<TH1F*>((*this)[indices]->Clone(name.data()));

    sum->Reset("MICES");
    for (int64_t i = 0; i < _shape[axis]; ++i) {
        indices[axis] = i;
        sum->Add((*this)[indices]);
    }

    return sum;
}

void history::apply(std::function<void(TH1*)> f) {
    for (auto& hist : histograms) { f(hist); }
}

void history::apply(std::function<void(TH1*, int64_t)> f) {
    for (int64_t i = 0; i < _size; ++i) { f(histograms[i], i); }
}

void history::save(std::string const& prefix) const {
    for (auto const& hist : histograms) {
        auto name = prefix + "_"s + hist->GetName();
        hist->Write(name.data(), TObject::kOverwrite);
    }

    auto shape_desc = ""s;
    for (auto const& s : _shape)
        shape_desc += "_"s + std::to_string(s);

    auto title = prefix + "_"s + _tag;
    auto label = new TNamed(title.data(), shape_desc.data());
    label->Write("", TObject::kOverwrite);
}

void history::rename(std::string const& prefix) {
    rename("", prefix);
}

void history::rename(std::string const& replace, std::string const& prefix) {
    auto original = _tag;

    _tag = prefix + "_"s + (replace.empty() ? _tag : replace);
    for (auto const& hist : histograms) {
        std::string name = hist->GetName();
        auto pos = name.find(original);
        name.replace(pos, original.length(), _tag);
        hist->SetName(name.data());
    }
}

bool history::compatible(history const& other) const {
    if (_dims < other._dims) { return false; }

    std::vector<int64_t> common = _shape;
    common.resize(other._dims);
    if (common != other._shape) { return false; }

    return true;
}

std::string history::stub(std::vector<int64_t> const& indices) const {
    auto add = [](std::string base, int64_t index) {
        return std::move(base) + "_"s + std::to_string(index); };

    return std::accumulate(std::begin(indices), std::end(indices), ""s, add);
}

history* history::_sum(int64_t axis) const {
    std::vector<int64_t> output = _shape;
    output.erase(std::next(std::begin(output), axis));

    auto result = new history(_tag + "_sum"s + std::to_string(axis),
                              _ordinate, output);

    for (int64_t i = 0; i < result->size(); ++i) {
        auto indices = result->indices_for(i);
        indices.insert(std::next(std::begin(indices), axis), 0);
        (*result)[i] = this->sum(indices, axis);
    }

    return result;
}

history* history::shrink(std::string const& tag,
        std::vector<int64_t> const& shape,
        std::vector<int64_t> const& offset) const {
    auto result = new history(*this, tag);

    auto pos = std::begin(result->histograms);
    for (int64_t i = 0; i < _size; ++i, ++pos) {
        auto const& indices = indices_for(i);
        for (int64_t j = 0; j < _dims; ++j) {
            if (indices[j] < offset[j] || indices[j] >= shape[j] + offset[j]) {
                pos = result->histograms.erase(pos);
                --pos;
                break;
            }
        }
    }

    result->_shape = shape;
    result->_size = std::accumulate(std::begin(shape), std::end(shape), 1,
                                    std::multiplies<int64_t>());

    result->apply([&](TH1* h, int64_t i) {
        h->SetName((result->_tag + stub(result->indices_for(i))).data()); });

    return result;
}

void history::allocate_histograms() {
    histograms = std::vector<TH1F*>(_size, nullptr);
    for (int64_t i = 0; i < _size; ++i) {
        histograms[i] = bins->book<TH1F>(_tag + stub(indices_for(i)),
            ";"s + bins->abscissa() + ";"s + _ordinate);
    }
}
