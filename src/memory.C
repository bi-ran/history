#include "../include/memory.h"

memory::memory(memory const& other, std::string const& prefix)
    : history(other, prefix),
      intervals(other.intervals) {
}

memory::memory(history&& other, std::shared_ptr<multival> const& intervals)
    : history(std::move(other)),
      intervals(intervals) {
}
