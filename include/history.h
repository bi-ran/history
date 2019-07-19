#ifndef HISTORY_H
#define HISTORY_H

#include "interval.h"
#include "multival.h"

#include "TFile.h"
#include "TH1.h"
#include "TObject.h"

#include <array>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using x = std::initializer_list<int64_t>;
using v = std::initializer_list<double>;

class history {
  public:
    history(std::string const& tag, std::string const& ordinate,
            std::shared_ptr<interval> const& bins,
            std::shared_ptr<multival> const& intervals)
            : _tag(tag),
              _ordinate(ordinate),
              _dims(intervals->dims()),
              _size(intervals->size()),
              _shape(intervals->shape()),
              bins(bins),
              intervals(intervals) {
        allocate_histograms();
    }

    template <template <typename...> class T>
    history(std::string const& tag, std::string const& ordinate,
            std::shared_ptr<interval> const& bins, T<int64_t> const& shape)
            : _tag(tag),
              _ordinate(ordinate),
              _dims(shape.size()),
              _size(std::accumulate(std::begin(shape), std::end(shape), 1,
                                    std::multiplies<int64_t>())),
              _shape(std::vector<int64_t>(std::begin(shape), std::end(shape))),
              bins(bins) {
        allocate_histograms();
    }

    template <typename... T>
    history(std::string const& tag, std::string const& ordinate,
            std::shared_ptr<interval> const& bins, T const&... dimensions)
            : _tag(tag),
              _ordinate(ordinate),
              _dims(sizeof...(T)),
              _size(size_of(dimensions...)),
              bins(bins) {
        auto shape = std::array<int64_t, sizeof...(T)>(dimensions...);
        _shape = std::vector<int64_t>(std::begin(shape), std::end(shape));

        allocate_histograms();
    }

    template <template <typename...> class T>
    history(std::string const& tag,
            std::string const& ordinate,
            T<float> const& edges,
            std::shared_ptr<multival> const& intervals)
        : history(tag, ordinate, std::make_shared<interval>(edges),
                  intervals) {
    }

    template <template <typename...> class T, template <typename...> class U>
    history(std::string const& tag,
            std::string const& ordinate,
            T<float> const& edges,
            U<int64_t> const& shape)
        : history(tag, ordinate, std::make_shared<interval>(edges), shape) {
    }

    template <template <typename...> class T, typename... U>
    history(std::string const& tag,
            std::string const& ordinate,
            T<float> const& edges,
            U const&... dimensions)
        : history(tag, ordinate, std::make_shared<interval>(edges),
                  dimensions...) {
    }

    template <template <typename...> class T>
    history(std::string const& tag,
            std::string const& ordinate,
            std::string const& abscissa,
            T<float> const& edges,
            std::shared_ptr<multival> const& intervals)
        : history(tag, ordinate, std::make_shared<interval>(edges, abscissa),
                  intervals) {
    }

    template <template <typename...> class T, template <typename...> class U>
    history(std::string const& tag, std::string const& ordinate,
            std::string const& abscissa, T<float> const& edges,
            U<int64_t> const& shape)
        : history(tag, ordinate, std::make_shared<interval>(edges, abscissa),
                  shape) {
    }

    template <template <typename...> class T, typename... U>
    history(std::string const& tag, std::string const& ordinate,
            std::string const& abscissa, T<float> const& edges,
            U const&... dimensions)
        : history(tag, ordinate, std::make_shared<interval>(edges, abscissa),
                  dimensions...) {
    }

    history(TFile* f, std::string const& tag);

    history(history const&, std::string const& prefix);

    history(history const&) = delete;
    history& operator=(history const&) = delete;
    ~history() = default;

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
        return intervals->index_for(values); }

    std::vector<int64_t> indices_for(int64_t index) const;

    void add(history const& other, double c1);

    void operator+=(history const& other);
    void operator-=(history const& other);

    void scale(double c1);

    void operator*=(double c1);
    void operator/=(double c1);

    void multiply(history const& other);
    void divide(history const& other);

    void operator*=(history const& other);
    void operator/=(history const& other);

    /* scale histograms integrated along axes. assume self, other have equal
     * shapes after integrating out axes */
    template <template <typename...> class T>
    void multiply(history const& other, T<int64_t> axes) {
        for (int64_t j = 0; j < other.size(); ++j) {
            auto count = other[j]->GetBinContent(1);
            auto indices = other.indices_for(j);
            for (auto const& axis : axes)
                indices.insert(std::next(std::begin(indices), axis), 0);

            std::function<void(std::vector<int64_t> const&)> scaler =
                    [&](std::vector<int64_t> const& indices) {
                (*this)[indices]->Scale(count); };

            permute(scaler, indices, _shape, axes);
        }
    }

    template <template <typename...> class T>
    void divide(history const& other, T<int64_t> axes) {
        for (int64_t j = 0; j < other.size(); ++j) {
            auto count = other[j]->GetBinContent(1);
            auto scale = count != 0 ? 1. / count : 0;
            auto indices = other.indices_for(j);
            for (auto const& axis : axes)
                indices.insert(std::next(std::begin(indices), axis), 0);

            std::function<void(std::vector<int64_t> const&)> scaler =
                    [&](std::vector<int64_t> const& indices) {
                (*this)[indices]->Scale(scale); };

            permute(scaler, indices, _shape, axes);
        }
    }

    void multiply(TH1* const other);
    void divide(TH1* const other);

    TH1F*& operator[](int64_t index);
    TH1F* const& operator[](int64_t index) const;

    template <template <typename...> class T, typename U>
    TH1F*& operator[](T<U> const& indices) {
        return histograms[index_for(indices)]; }

    template <template <typename...> class T, typename U>
    TH1F* const& operator[](T<U> const& indices) const {
        return histograms[index_for(indices)]; }

    TH1F* sum(std::vector<int64_t> const& indices, int64_t axis) const;

    TH1F* sum(std::vector<int64_t> const& indices, int64_t axis,
              int64_t start, int64_t end) const;

    std::unique_ptr<history> sum(int64_t axis) const {
        return sum_impl(axis); }

    template <typename... T>
    std::unique_ptr<history> sum(int64_t axis, T... axes) const {
        return sum_impl(axis)->sum(axes...); }

    template <typename T, typename... U>
    T operator()(int64_t index, T (TH1::* function)(U...), U... args) {
        return forward(index, function, args...); }

    template <typename T, typename... U>
    T operator()(int64_t index, T (TH1::* function)(U...) const,
                 U... args) const {
        return forward(index, function, args...); }

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (TH1::* function)(W...),
                 W... args) {
        return forward(index_for(indices), function, args...); }

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (TH1::* function)(W...) const,
                 W... args) {
        return forward(index_for(indices), function, args...); }

    void apply(std::function<void(TH1*)> f);

    void save(std::string const& prefix) const;

    int64_t const& dims() const { return _dims; }
    int64_t const& size() const { return _size; }
    std::vector<int64_t> const& shape() const { return _shape; }

  private:
    void _multiply(history const& other);
    void _divide(history const& other);

    template <template <typename...> class T, typename U>
    void permute(std::function<void(T<U> const&)>& lambda,
                 std::vector<int64_t>& indices,
                 std::vector<int64_t> const& shape,
                 std::vector<int64_t> const& axes,
                 int64_t index = 0) {
        if (index == static_cast<int64_t>(axes.size())) {
            lambda(indices); return; }

        int64_t axis = axes[index];
        for (int64_t i = 0; i < shape[axis]; ++i) {
            indices[axis] = i;
            permute(lambda, indices, shape, axes, index + 1);
        }
    }

    TH1F* sum_impl(std::string const& name, std::vector<int64_t> indices,
                   int64_t axis, int64_t start, int64_t end) const;

    std::unique_ptr<history> sum_impl(int64_t axis) const;

    template <typename T, typename U, typename... V>
    T forward(int64_t index, T (U::* function)(V...), V... args) {
        return ((*histograms[index]).*function)(std::forward<V>(args)...); }

    template <typename T, typename U, typename... V>
    T forward(int64_t index, T (U::* function)(V...) const, V... args) const {
        return ((*histograms[index]).*function)(std::forward<V>(args)...); }

    void allocate_histograms();

    template <typename... T>
    constexpr int64_t size_of(T const&... dimensions) const {
        int64_t size = 1;
        for (auto dim : { dimensions... })
            size = size * dim;
        return size;
    }

    std::string _tag;
    std::string _ordinate;

    int64_t _dims;
    int64_t _size;
    std::vector<int64_t> _shape;

    std::shared_ptr<interval> bins;
    std::shared_ptr<multival> intervals;
    std::vector<TH1F*> histograms;
};

#endif /* HISTORY_H */
