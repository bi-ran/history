#ifndef HISTORY_H
#define HISTORY_H

#include "TFile.h"
#include "TNamed.h"
#include "TObject.h"

#include <array>
#include <functional>
#include <iterator>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using x = std::initializer_list<int64_t> const;
using v = std::initializer_list<double> const;

template <typename H>
class history {
  public:
    history(std::string const& tag, std::string const& label,
            std::vector<int64_t> const& shape)
        : _tag(tag),
          _label(label),
          _dims(shape.size()),
          _size(std::accumulate(std::begin(shape), std::end(shape), 1,
                                std::multiplies<int64_t>())),
          _shape(shape),
          objects(std::vector<H*>(_size, nullptr)) {
    }

    template <template <typename...> class T>
    history(std::string const& tag, std::string const& label,
            auto factory, T<int64_t> const& shape)
            : _tag(tag),
              _label(label),
              _dims(shape.size()),
              _size(std::accumulate(std::begin(shape), std::end(shape), 1,
                                    std::multiplies<int64_t>())),
              _shape(std::vector<int64_t>(std::begin(shape), std::end(shape))),
              _factory(factory) {
        allocate_objects();
    }

    template <typename... T>
    history(std::string const& tag, std::string const& label,
            auto factory, T const&... dimensions)
            : _tag(tag),
              _label(label),
              _dims(sizeof...(T)),
              _size(size_of(dimensions...)),
              _factory(factory) {
        auto shape = std::array<int64_t, sizeof...(T)>({dimensions...});
        _shape = std::vector<int64_t>(std::begin(shape), std::end(shape));

        allocate_objects();
    }

    history(TFile* f, std::string const& tag)
            : _tag(tag) {
        std::string desc = ((TNamed*)f->Get(tag.data()))->GetTitle();
        while (!desc.empty()) {
            desc.erase(0, 1);

            auto pos = desc.find("_");
            auto token = desc.substr(0, pos);
            _shape.push_back(std::atoi(token.data()));

            desc.erase(0, pos);
        }

        _dims = _shape.size();
        _size = std::accumulate(std::begin(_shape), std::end(_shape), 1,
                                std::multiplies<int64_t>());

        objects = std::vector<H*>(_size, nullptr);
        for (int64_t i = 0; i < _size; ++i) {
            auto name = _tag + stub(i);
            objects[i] = (H*)f->Get(name.data());
            objects[i]->SetName(name.data());
        }
    }

    template <typename... T>
    history(TFile* f, std::string const& tag, T const&... args)
            : history(f, tag) {
        rename(args...);
    }

    history(history const& other, std::string const& prefix)
            : _tag(prefix + "_" + other._tag),
              _label(other._label),
              _dims(other._dims),
              _size(other._size),
              _shape(other._shape),
              _factory(other._factory) {
        for (auto const& obj : other.objects)
            objects.push_back((H*)obj->Clone());

        rename();
    }

    history(history const& other, std::string const& old,
            std::string const& tag)
            : history(other, "") {
        rename(std::string("_") + old, tag);
    }

    history(history const&) = delete;
    history& operator=(history const&) = delete;
    history(history&&) = default;
    history& operator=(history&&) = default;
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

    std::vector<int64_t> indices_for(int64_t index) const {
        std::vector<int64_t> indices(_dims);
        for (int64_t i = 0; i < _dims; ++i) {
            indices[i] = index % _shape[i];
            index = index / _shape[i];
        }

        return indices;
    }

    void add(history const& other, double c1) {
        for (int64_t i = 0; i < _size; ++i)
            objects[i]->Add(other[i], c1);
    }

    void operator+=(history const& other) { add(other, 1); }
    void operator-=(history const& other) { add(other, -1); }

    void scale(double c1) {
        for (auto const& obj : objects)
            obj->Scale(c1);
    }

    void operator*=(double c1) { scale(c1); }
    void operator/=(double c1) { scale(1. / c1); }

    /* scale histograms by counts from other. assume self, other have identical
     * shapes for overlapping dimensions */
    void multiply(history const& other) {
        /* warning: fails silently! */
        if (!compatible(other)) { return; }

        std::vector<int64_t> axes(_dims - other._dims);
        std::iota(std::begin(axes), std::end(axes), other._dims);

        return multiply(other, axes);
    }

    void divide(history const& other) {
        /* warning: fails silently! */
        if (!compatible(other)) { return; }

        std::vector<int64_t> axes(_dims - other._dims);
        std::iota(std::begin(axes), std::end(axes), other._dims);

        return divide(other, axes);
    }

    void operator*=(history const& other) { multiply(other); }
    void operator/=(history const& other) { divide(other); }

    /* scale histograms integrated along axes. assume self, other have equal
     * shapes after integrating out axes */
    template <template <typename...> class T>
    void multiply(history const& other, T<int64_t> axes) {
        scale(other, axes, [](float content) -> float {
            return content; });
    }

    template <template <typename...> class T>
    void divide(history const& other, T<int64_t> axes) {
        scale(other, axes, [](float content) -> float {
            return content != 0 ? 1. / content : 0; });
    }

    H*& operator[](int64_t index) { return objects[index]; }
    H* const& operator[](int64_t index) const { return objects[index]; }

    template <template <typename...> class T, typename U>
    H*& operator[](T<U> const& indices) {
        return objects[index_for(indices)]; }

    template <template <typename...> class T, typename U>
    H* const& operator[](T<U> const& indices) const {
        return objects[index_for(indices)]; }

    H* sum(std::vector<int64_t> indices, int64_t axis) const {
        std::vector<int64_t> output = indices;
        output.erase(std::next(std::begin(output), axis));

        auto name = _tag + "_sum" + std::to_string(axis) + stub(output);
        auto sum = static_cast<H*>((*this)[indices]->Clone(name.data()));

        sum->Reset("MICES");
        for (int64_t i = 0; i < _shape[axis]; ++i) {
            indices[axis] = i;
            sum->Add((*this)[indices]);
        }

        return sum;
    }

    history* sum(int64_t axis) const {
        std::vector<int64_t> output = _shape;
        output.erase(std::next(std::begin(output), axis));

        auto result = new history(_tag + "_sum" + std::to_string(axis),
                                  _label, output);

        for (int64_t i = 0; i < result->size(); ++i) {
            auto indices = result->indices_for(i);
            indices.insert(std::next(std::begin(indices), axis), 0);
            (*result)[i] = this->sum(indices, axis);
        }

        return result;
    }

    template <typename... T>
    history* sum(int64_t axis, T... axes) const {
        return sum(axis)->sum(axes...); }

    history* extend(std::string const& prefix, int64_t axis, int64_t size) {
        auto shape = std::vector<int64_t>(_shape);
        shape.insert(std::next(std::begin(shape), axis), size);

        auto result = new history(prefix + "_" + _tag, _label, shape);

        for (int64_t i = 0; i < _size; ++i) {
            auto indices = indices_for(i);
            indices.insert(std::next(std::begin(indices), axis), 0);

            for (int64_t j = 0; j < size; ++j) {
                indices[axis] = j;
                (*result)[indices] = (H*)objects[i]->Clone();
            }
        }

        result->rename();

        return result;
    }

    history* shrink(std::string const& tag,
                    std::vector<int64_t> const& shape,
                    std::vector<int64_t> const& offset) const {
        auto result = new history(*this, tag);

        std::vector<int64_t> out;
        std::function<void(std::vector<int64_t> const&)> recorder =
                [&](std::vector<int64_t> const& indices) {
            for (int64_t j = 0; j < _dims; ++j) {
                if (indices[j] >= shape[j] + offset[j]
                        || indices[j] < offset[j])
                    out.push_back(index_for(indices));
            }
        };

        std::vector<int64_t> zero(_dims, 0);
        std::vector<int64_t> axes(_dims);
        std::iota(std::begin(axes), std::end(axes), 0);

        permute(recorder, zero, _shape, axes, 0);

        auto ref = std::begin(result->objects);
        for (auto it = std::rbegin(out); it != std::rend(out); ++it)
            result->objects.erase(std::next(ref, *it));

        result->_shape = shape;
        result->_size = std::accumulate(std::begin(shape), std::end(shape), 1,
                                        std::multiplies<int64_t>());

        result->rename();

        return result;
    }

    template <typename T, typename... U>
    T operator()(int64_t index, T (H::* fn)(U...), U... args) {
        return forward(index, fn, args...); }

    template <typename T, typename... U>
    T operator()(int64_t index, T (H::* fn)(U...) const, U... args) const {
        return forward(index, fn, args...); }

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (H::* fn)(W...), W... args) {
        return forward(index_for(indices), fn, args...); }

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (H::* fn)(W...) const,
                 W... args) const {
        return forward(index_for(indices), fn, args...); }

    void apply(std::function<void(H*)> f) {
        for (auto& obj : objects) { f(obj); } }

    void apply(std::function<void(H*, int64_t)> f) {
        for (int64_t i = 0; i < _size; ++i) { f(objects[i], i); } }

    void save(std::string const& prefix) const {
        auto full = prefix.empty() ? "" : prefix + "_";
        for (auto const& obj : objects) {
            auto name = full + obj->GetName();
            obj->Write(name.data(), TObject::kOverwrite);
        }

        auto label = new TNamed((full + _tag).data(), stub(_shape).data());
        label->Write("", TObject::kOverwrite);
    }

    void save() const { save(""); }

    void rename(std::string const& old, std::string const& tag) {
        _tag.replace(_tag.find(old), old.length(), tag); rename(); }

    void rename(std::string const& tag) {
        _tag = tag; rename(); }

    void rename() {
        for (int64_t i = 0; i < _size; ++i)
            objects[i]->SetName((_tag + stub(i)).data());
    }

    int64_t const& dims() const { return _dims; }
    int64_t const& size() const { return _size; }
    std::vector<int64_t> const& shape() const { return _shape; }

  protected:
    bool compatible(history const& other) const {
        if (_dims < other._dims) { return false; }

        std::vector<int64_t> common = _shape;
        common.resize(other._dims);
        if (common != other._shape) { return false; }

        return true;
    }

    template <template <typename...> class T, typename U>
    void permute(std::function<void(T<U> const&)>& lambda,
                 std::vector<int64_t>& indices,
                 std::vector<int64_t> const& shape,
                 std::vector<int64_t> const& axes,
                 int64_t index) const {
        if (index == static_cast<int64_t>(axes.size())) {
            lambda(indices); return; }

        int64_t axis = axes[index];
        for (int64_t i = 0; i < shape[axis]; ++i) {
            indices[axis] = i;
            permute(lambda, indices, shape, axes, index + 1);
        }
    }

    template <template <typename...> class T>
    void scale(history const& other, T<int64_t> axes,
               std::function<float(float)> lambda) {
        for (int64_t j = 0; j < other.size(); ++j) {
            auto indices = other.indices_for(j);
            for (auto const& axis : axes)
                indices.insert(std::next(std::begin(indices), axis), 0);

            auto scale = lambda(other[j]->GetBinContent(1));
            std::function<void(std::vector<int64_t> const&)> scaler =
                    [&](std::vector<int64_t> const& indices) {
                (*this)[indices]->Scale(scale); };

            permute(scaler, indices, _shape, axes, 0);
        }
    }

    std::string stub(std::vector<int64_t> const& indices) const {
        return std::accumulate(std::begin(indices), std::end(indices),
            std::string(), [](std::string base, int64_t index) {
                return std::move(base) + "_" + std::to_string(index); });
    }

    std::string stub(int64_t index) const {
        return stub(indices_for(index)); }

    template <typename T, typename... U>
    T forward(int64_t index, T (H::* function)(U...), U... args) {
        return ((*objects[index]).*function)(std::forward<U>(args)...); }

    template <typename T, typename... U>
    T forward(int64_t index, T (H::* function)(U...) const, U... args) const {
        return ((*objects[index]).*function)(std::forward<U>(args)...); }

    void allocate_objects() {
        objects = std::vector<H*>(_size, nullptr);
        for (int64_t i = 0; i < _size; ++i)
            objects[i] = _factory(i, _tag + stub(i), _label);
    }

    template <typename... T>
    constexpr int64_t size_of(T const&... dimensions) const {
        auto dims = { dimensions... };
        return std::accumulate(std::begin(dims), std::end(dims), 1,
                               std::multiplies<int64_t>());
    }

    std::string _tag;
    std::string _label;

    int64_t _dims;
    int64_t _size;
    std::vector<int64_t> _shape;

    std::function<H*(int64_t, std::string const&, std::string const&)> _factory;
    std::vector<H*> objects;
};

#endif /* HISTORY_H */
