#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>
#include <vector>

class interval {
  public:
    interval(int64_t number, double min, double max,
             std::string const& abscissa);

    interval(int64_t number, double min, double max);

    template <template <typename...> class T>
    interval(T<float> const& edges, std::string const& abscissa)
        : _size(edges.size() - 1),
          _edges(std::begin(edges), std::end(edges)),
          _abscissa(abscissa) { }

    template <template <typename...> class T>
    interval(T<float> const& edges)
        : interval(edges, std::string()) { }

    template <typename... T>
    interval(T const&... edges, std::string const& abscissa)
        : _size(sizeof...(T) - 1),
          _edges({static_cast<double>(edges)...}),
          _abscissa(abscissa) { }

    template <typename... T>
    interval(T const&... edges)
        : _size(sizeof...(T) - 1),
          _edges({static_cast<double>(edges)...}) { }

    interval(interval const& other) = default;
    interval& operator=(interval const& other) = default;
    ~interval() = default;

    int64_t index_for(double value) const;

    double operator[](int64_t index) const { return _edges[index]; }

    int64_t const& size() const { return _size; }
    std::string const& abscissa() const { return _abscissa; }

    template <typename T>
    T* book(std::string const& name, std::string const& title);

  private:
    int64_t _size;
    std::vector<double> _edges;

    std::string const _abscissa;
};

#endif /* INTERVAL_H */
