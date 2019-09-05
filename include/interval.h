#ifndef INTERVAL_H
#define INTERVAL_H

#include <array>
#include <string>
#include <vector>

class interval {
  public:
    interval(std::string const& abscissa, int64_t number,
             double min, double max);

    interval(int64_t number, double min, double max);

    template <template <typename...> class T>
    interval(std::string const& abscissa, T<float> const& edges)
        : _abscissa(abscissa),
          _size(edges.size() - 1),
          _edges(std::begin(edges), std::end(edges)) { }

    template <template <typename...> class T>
    interval(T<float> const& edges)
        : interval(std::string(), edges) { }

    interval(interval const& other) = default;
    interval& operator=(interval const& other) = default;
    ~interval() = default;

    int64_t index_for(double value) const;

    double const* edges() const {
        return _edges.data(); }
    std::array<double, 2> edges(int64_t index) const {
        return { _edges[index], _edges[index + 1] }; }

    double operator[](int64_t index) const {
        return _edges[index]; }

    std::string const& abscissa() const { return _abscissa; }
    int64_t size() const { return _size; }

  private:
    std::string const _abscissa;

    int64_t const _size;
    std::vector<double> _edges;
};

#endif /* INTERVAL_H */
