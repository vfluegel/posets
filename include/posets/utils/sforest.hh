#pragma once

#include <cassert>
#include <iostream>
#include <ranges>
#include <stack>
#include <unordered_map>
#include <vector>

#include <posets/concepts.hh>

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

  struct st_hash {
    size_t operator()(const st_node& k) const {
      size_t res = std::hash<int>()(k.label);
      for (size_t i = 0; i < k.numchild; i++)
        res ^= std::hash<size_t>()(k.children[i]) << (i + 1);
      return res;
    }
  };

  struct st_equal {
    bool operator()(const st_node& lhs, const st_node& rhs) const {
      if (lhs.label != rhs.label or lhs.numchild != rhs.numchild) return false;
      for (size_t i = 0; i < lhs.numchild; i++)
        if (lhs.children[i] != rhs.children[i]) return false;
      return true;
    }
  };

  st_node **layers;
  size_t *layer_size;
  size_t *layer_nxt;
  size_t *child_buffer;
  size_t cbuffer_size;
  size_t cbuffer_nxt;

  std::unordered_map<st_node, size_t, st_hash, st_equal> *inverse;

  void init(int k, size_t dim) {
    this->k = k;
    this->dim = dim;
    layers = new st_node*[dim];
    inverse = new std::unordered_map<st_node, size_t, st_hash, st_equal>[dim];
    for (size_t i = 0; i < dim; i++) {
      layer_size[i] = INIT_LAYER_SIZE;
      layers[i] = new st_node[layer_size[i]];
      layer_nxt[i] = 0;
    }
    cbuffer_size = INIT_LAYER_SIZE * (k + 1);
    child_buffer = new size_t[cbuffer_size];
    cbuffer_nxt = 0;
  }

public:
  sforest() : layers{nullptr} {}

  sforest(int k, size_t dim) {
    this->init(k, dim);
  }

  ~sforest() {
    if (layers == nullptr) return;

    for (size_t i = 0; i < dim; i++)
      delete[] layers[i];
    delete[] layers;
    delete[] layer_size;
    delete[] layer_nxt;
    delete[] inverse;
    delete[] child_buffer;
  }

  std::vector<V> get_all() {
    // Declare a stack of tuples of stnode, child no., and layer
    // Traverse all roots at dimension 0
    for (size_t i = 0; i < layer_nxt[0]; i++) {
      // TODO
    }
  }

  template <std::ranges::input_range R>
  void add_vectors(R&& elements) {
    // TODO: Vanessa please help
  }
};

// FIXME: This has to be built recursively from the forest
template <Vector V>
inline std::ostream& operator<< (std::ostream& os, const sforest<V>& f) {
  for (auto&& el : f.get_all())
    os << el << std::endl;

  return os;
}

} // namespace posets::utils
