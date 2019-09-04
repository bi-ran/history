#ifndef MEMORY_H
#define MEMORY_H

#include "history.h"

#include "multival.h"

class memory : public history {
  public:
    memory(std::string const& tag, std::string const& ordinate,
           interval const* bins, multival const* intervals)
        : history(tag, ordinate, bins, intervals->shape()),
          intervals(intervals) {
    }

    template <template <typename...> class T>
    memory(std::string const& tag,
           std::string const& ordinate,
           T<float> const& edges,
           multival const* intervals)
        : history(tag, ordinate, new interval(edges), intervals->shape()),
          intervals(intervals) {
    }

    template <template <typename...> class T>
    memory(std::string const& tag,
           std::string const& ordinate,
           std::string const& abscissa,
           T<float> const& edges,
           multival const* intervals)
        : history(tag, ordinate, new interval(abscissa, edges),
                  intervals->shape()),
          intervals(intervals) {
    }

    memory(history&&, multival const* intervals);

    memory(memory&&) = delete;
    memory& operator=(memory&&) = delete;

    using history::index_for;

    template <template <typename...> class T, typename U>
    typename std::enable_if<std::is_floating_point<U>::value, int64_t>::type
    index_for(T<U> const& values) const {
        return intervals->index_for(values); }

    using history::operator[];

    template <template <typename...> class T, typename U>
    TH1F*& operator[](T<U> const& indices) {
        return histograms[index_for(indices)]; }

    template <template <typename...> class T, typename U>
    TH1F* const& operator[](T<U> const& indices) const {
        return histograms[index_for(indices)]; }

    using history::operator();

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (TH1::* fn)(W...), W... args) {
        return forward(index_for(indices), fn, args...); }

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (TH1::* fn)(W...) const,
                 W... args) const {
        return forward(index_for(indices), fn, args...); }

    memory(memory const&, std::string const& prefix);

    memory(memory const&) = delete;
    memory& operator=(memory const&) = delete;
    ~memory() = default;

  private:
    multival const* intervals;
};

#endif /* MEMORY_H */
