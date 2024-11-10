#pragma once

#include <algorithm>
#include <cassert>
#include <map>
#include <numeric>
#include <ranges>
#include <vector>

namespace posets::utils {

size_t INIT_LAYER_SIZE = 100;

// Forward definition for the operator<<
template <Vector>
class sforest;

template <Vector V>
std::ostream& operator<< (std::ostream& os, const utils::sforest<V>& f);

template <Vector V> class sforest {
private:
  int k;
  size_t dim;

  struct st_node {
    int label;
    size_t *children;
    size_t numchild;
  };

  st_node **layers;
  size_t *child_buffer;
  size_t cbuffer_size;
  size_t cbuffer_nxt;

  std::map<std::vector<size_t>, size_t> cache;

public:
  template <std::ranges::input_range R>
  sharingforest(int k, size_t dim)
      : k{k}, dim{dim} {
    layers = new st_node*[dim];
    for (size_t i = 0; i < dim; i++)
      layers[i] = new st_node[INIT_LAYER_SIZE];
    cbuffer_size = INIT_LAYER_SIZE * (k + 1);
    child_buffer = new size_t[cbuffer_size];
    cbuffer_nxt = 0;
  }

  ~sharingforest() {
    for (size_t i = 0; i < dim; i++)
      delete[] layers[i];
    delete[] layers;
    delete[] child_buffer;
  }
};

// FIXME: This has to be built recursively from the forest
template <Vector V>
inline std::ostream& operator<< (std::ostream& os, const sharingforest<V>& f) {
  for (auto&& el : f.vector_set)
    os << el << std::endl;

  return os;
}

} // namespace posets::utils
