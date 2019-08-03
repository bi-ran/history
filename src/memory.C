#include "../include/memory.h"

memory::memory(memory const& other, std::string const& prefix)
    : history(other, prefix),
      intervals(other.intervals) {
}
