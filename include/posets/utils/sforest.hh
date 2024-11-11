#pragma once

#include <cassert>
#include <iostream>
#include <ranges>
#include <stack>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <posets/concepts.hh>

namespace posets::utils {

size_t INIT_LAYER_SIZE = 100;

// Forward definition for the operator<<
template <Vector> class sforest;

template <Vector V>
std::ostream &operator<<(std::ostream &os, const utils::sforest<V> &f);

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
    size_t operator()(const st_node &k) const {
      size_t res = std::hash<int>()(k.label);
      for (size_t i = 0; i < k.numchild; i++)
        res ^= std::hash<size_t>()(k.children[i]) << (i + 1);
      return res;
    }
  };

  struct st_equal {
    bool operator()(const st_node &lhs, const st_node &rhs) const {
      if (lhs.label != rhs.label or lhs.numchild != rhs.numchild)
        return false;
      for (size_t i = 0; i < lhs.numchild; i++)
        if (lhs.children[i] != rhs.children[i])
          return false;
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
    layers = new st_node *[dim];
    layer_size = new size_t[dim];
    layer_nxt = new size_t[dim];
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
    assert(dim >= 2);
    this->init(k, dim);
  }

  ~sforest() {
    if (layers == nullptr)
      return;

    for (size_t i = 0; i < dim; i++)
      delete[] layers[i];
    delete[] layers;
    delete[] layer_size;
    delete[] layer_nxt;
    delete[] inverse;
    delete[] child_buffer;
  }

  std::vector<V> get_all() {
    // Stack with tuples (layer, node id, child id)
    std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

    // Add all roots at dimension 0
    for (size_t i = 0; i < layer_nxt[0]; i++)
      to_visit.push({0, i, 0});

    std::vector<V> res;
    std::vector<int> temp;
    while (to_visit.size() > 0) {
      const auto [lay, node, child] = to_visit.top();
      to_visit.pop();
      const auto parent = layers[lay][node];
      // first time we see parent? get its label
      if (child == 0)
        temp.push_back(parent.label);

      // base case: reached the bottom layer
      if (lay == this->dim - 2) {
        assert(child == 0);
        for (size_t i = 0; i < parent.numchild; i++) {
          auto bottom_node = layers[lay + 1][parent.children[i]];
          temp.push_back(bottom_node.label);
          std::vector<int> cpy{temp};
          res.push_back(V(std::move(cpy)));
          temp.pop_back();
        }
        temp.pop_back(); // done with this parent
      } else {           // recursive case
        // Either we're done with this node and we just mark it as visited or
        // we need to keep it and we add it's next son
        if (child < parent.numchild) {
          to_visit.push({lay, node, child + 1});
          to_visit.push({lay + 1, parent.children[child], 0});
        } else {
          temp.pop_back(); // done with this parent
        }
      }
    }
  }

  bool cover_vector(std::vector<size_t> roots, V covered) {
    // Stack with tuples (layer, node id, child id)
    std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

    // Add all roots at dimension 0 such that their labels cover the first
    // component of the given vector
    for (const auto i : roots) {
      assert(i < layer_nxt[0]);
      if (covered[0] <= layers[0][i].label)
        to_visit.push({0, i, 0});
    }

    while (to_visit.size() > 0) {
      const auto [lay, node, child] = to_visit.top();
      to_visit.pop();
      const auto parent = layers[lay][node];

      // base case: reached the bottom layer
      if (lay == this->dim - 2) {
        assert(child == 0);
        // it is sufficient to check the last child
        size_t i = parent.numchild - 1;
        auto bottom_node = layers[lay + 1][parent.children[i]];
        if (covered[lay + 1] <= bottom_node.label)
          return true;
      } else { // recursive case
        // Either we're done with this node and we just mark it as visited or
        // we need to keep it and we add it's next son
        size_t c = child;
        if (c == 0) {
          size_t i = parent.numchild - 1;
          auto child_node = layers[lay + 1][parent.children[i]];
          // find the first index where the order holds
          // FIXME: this could be a binary search
          if (covered[lay] > child_node.label)
            continue;
          do {
            child_node = layers[lay + 1][parent.children[c]];
            c++;
          } while (covered[lay] > child_node.label);
          c--;
        }
        if (c < parent.numchild) {
          to_visit.push({lay, node, c + 1});
          to_visit.push({lay + 1, parent.children[c], 0});
        }
      }
    }
    return false;
  }

  template <std::ranges::input_range R>
  std::vector<size_t> add_vectors(R &&elements) {
    // TODO: Vanessa please help, this should return the indices of the root
    // nodes (in the first layer) of the set of vectors added
    std::vector<size_t> res;
    return res;
  }
};

// FIXME: This has to be built recursively from the forest
template <Vector V>
inline std::ostream &operator<<(std::ostream &os, const sforest<V> &f) {
  for (auto &&el : f.get_all())
    os << el << std::endl;

  return os;
}

} // namespace posets::utils
