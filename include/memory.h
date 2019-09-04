#ifndef MEMORY_H
#define MEMORY_H

#include "history.h"
#include "multival.h"

template <typename H>
class memory : public history<H> {
  public:
    memory(std::string const& tag, std::string const& ordinate,
           interval const* bins, multival const* intervals)
        : history<H>(tag, ordinate, bins, intervals->shape()),
          intervals(intervals) {
    }

    template <template <typename...> class T>
    memory(std::string const& tag,
           std::string const& ordinate,
           T<float> const& edges,
           multival const* intervals)
        : history<H>(tag, ordinate, new interval(edges), intervals->shape()),
          intervals(intervals) {
    }

    template <template <typename...> class T>
    memory(std::string const& tag,
           std::string const& ordinate,
           std::string const& abscissa,
           T<float> const& edges,
           multival const* intervals)
        : history<H>(tag, ordinate, new interval(abscissa, edges),
                  intervals->shape()),
          intervals(intervals) {
    }

    memory(memory const& other, std::string const& prefix)
        : history<H>(other, prefix),
          intervals(other.intervals) {
    }

    memory(history<H>&& other, multival const* intervals)
        : history<H>(std::move(other)),
          intervals(intervals) {
    }

    memory(memory const&) = delete;
    memory& operator=(memory const&) = delete;
    memory(memory&&) = delete;
    memory& operator=(memory&&) = delete;
    ~memory() = default;

    using history<H>::index_for;

    template <template <typename...> class T, typename U>
    typename std::enable_if<std::is_floating_point<U>::value, int64_t>::type
    index_for(T<U> const& values) const {
        return intervals->index_for(values); }

    using history<H>::operator[];

    template <template <typename...> class T, typename U>
    H*& operator[](T<U> const& indices) {
        return this->histograms[index_for(indices)]; }

    template <template <typename...> class T, typename U>
    H* const& operator[](T<U> const& indices) const {
        return this->histograms[index_for(indices)]; }

    using history<H>::operator();

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (H::* fn)(W...), W... args) {
        return forward(index_for(indices), fn, args...); }

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (H::* fn)(W...) const,
                 W... args) const {
        return forward(index_for(indices), fn, args...); }

  private:
    multival const* intervals;
};

#endif /* MEMORY_H */
