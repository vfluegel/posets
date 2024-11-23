#pragma once

#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <ranges>
#include <stack>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <posets/concepts.hh>

namespace posets::utils {

// We will be (micro)managing (C-style) dynamic memory for the set of nodes we
// keep in each layer. The initial number of nodes allowed is determined by
// the number below (multiplied by the number of layers = the dimension + 1).
// When the buffer is full, we double its capacity and copy everything to the
// new block of reserved memory.
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
    size_t numchild;
    size_t cbuffer_offset;
  };

  // NOTE: We will be keeping a unique table of st_nodes using an unordered
  // map, for this we need a hash function and an equivalence operation.
  // Since nodes keep offsets instead of pointers, we need both the hash and
  // equivalence operations to be able to access the unique table. Hence, we
  // keep a reference to the table in them.
  struct st_hash {
    sforest *f;
    st_hash(sforest *that) : f{that} {}

    size_t operator()(const st_node &k) const {
      size_t res = std::hash<int>()(k.label);
      size_t *children = f->child_buffer + k.cbuffer_offset;
      for (size_t i = 0; i < k.numchild; i++)
        res ^= std::hash<size_t>()(children[i]) << (i + 1);
      return res;
    }
  };

  struct st_equal {
    sforest *f;
    st_equal(sforest *that) : f{that} {}

    bool operator()(const st_node &lhs, const st_node &rhs) const {
      if (lhs.label != rhs.label or lhs.numchild != rhs.numchild)
        return false;
      size_t *lhs_children = f->child_buffer + lhs.cbuffer_offset;
      size_t *rhs_children = f->child_buffer + rhs.cbuffer_offset;
      for (size_t i = 0; i < lhs.numchild; i++)
        if (lhs_children[i] != rhs_children[i])
          return false;
      return true;
    }
  };

  std::vector<std::vector<st_node>> layers;
  size_t *child_buffer;
  size_t cbuffer_size;
  size_t cbuffer_nxt;

  std::vector<std::unordered_map<st_node, size_t, st_hash, st_equal>> inverse;

  void init(int k, size_t dim) {
    this->k = k;
    this->dim = dim;
    layers.resize(dim + 1);

    for (size_t i = 0; i < dim + 1; i++)
      inverse.emplace_back(INIT_LAYER_SIZE, st_hash(this), st_equal(this));

    cbuffer_size = INIT_LAYER_SIZE * (k + 1);
    child_buffer = new size_t[cbuffer_size];
    cbuffer_nxt = 0;
  }

  std::optional<size_t> hasSon(st_node &node, size_t childLayer, int val) {
    if (node.numchild == 0)
      return std::nullopt;

    size_t left = 0;
    size_t right = node.numchild - 1;
    size_t *children = child_buffer + node.cbuffer_offset;

    while (left <= right) {
      size_t mid = left + (right - left) / 2;
      assert(mid < node.numchild);
      int midVal = layers[childLayer][children[mid]].label;

      if (midVal == val) {
        return children[mid];
      } else if (midVal < val) {
        left = mid + 1;
      } else {
        right = mid - 1;
      }
    }

    return std::nullopt;
  }

  void addSon(st_node &node, size_t sonLayer, size_t son) {
    int last = node.numchild - 1;
    st_node &sonNode = layers[sonLayer][son];
    size_t *children = child_buffer + node.cbuffer_offset;

    // the new node is larger than all existing ones
    assert(last == -1 or
           layers[sonLayer][children[last]].label < sonNode.label);

    // Insert the new value
    children[last + 1] = son;
    node.numchild++;
  }

  size_t addNode(st_node &node, size_t destinationLayer) {
    auto existingNode = inverse[destinationLayer].find(node);
    if (existingNode == inverse[destinationLayer].end()) {
      size_t newID = layers[destinationLayer].size();
      inverse[destinationLayer][node] = newID;
      layers[destinationLayer].push_back(node);
      return newID;
    } else {
      return existingNode->second;
    }
  }

  size_t addChildren() {
    // Double the child buffer if it is full
    if (cbuffer_nxt + k >= cbuffer_size) {
      size_t *newBuffer = new size_t[cbuffer_size * 2];
      for (size_t i = 0; i < cbuffer_size; i++) {
        newBuffer[i] = child_buffer[i];
      }
      cbuffer_size *= 2;
      delete[] child_buffer;
      child_buffer = newBuffer;
    }
    size_t res = cbuffer_nxt;
    cbuffer_nxt += k;
    return res;
  }

  size_t build_node(std::vector<size_t> &vecs, size_t currentLayer,
                    const auto &elementVec) {
    assert(vecs.size() > 0);
    // If currentLayer is 0, we set the label to the dummy value -1 for the root
    // Else all nodes should have the same value at index currentLayer - 1, so
    // we just use the first
    int label{currentLayer == 0
                  ? -1
                  : static_cast<int>(elementVec[vecs[0]][currentLayer - 1])};
    st_node newNode{label, 0};
    // We have not reached the last layer - so add children
    if (currentLayer < this->dim) {
      newNode.cbuffer_offset = addChildren();
      // TODO: Try to do constact-access bucketing based on the value of k.
      // Probably won't pay off unless the set of vectors we are adding is
      // dense in most components.
      std::map<typename V::value_type, std::vector<size_t>>
          newPartition{}; // Partition and order the future children
      for (auto const &vec : vecs) {
        newPartition[elementVec[vec][currentLayer]].push_back(vec);
      }

      for (auto &[n, children] : newPartition) {
        // Build a new son for each individual value at currentLayer + 1
        size_t newSon = build_node(children, currentLayer + 1, elementVec);
        addSon(newNode, currentLayer + 1, newSon);
      }
    }
    return addNode(newNode, currentLayer);
  }

public:
  sforest() {}

  sforest(int k, size_t dim) {
    assert(dim >= 2);
    this->init(k, dim);
  }

  ~sforest() {
    if (layers.empty())
      return;

    delete[] child_buffer;
  }

  std::vector<V> get_all() const {
    // Stack with tuples (layer, node id, child id)
    std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

    // Add all roots at dimension 0
    for (size_t i = 0; i < layers[1].size(); i++)
      to_visit.push({1, i, 0});

    std::vector<V> res;
    std::vector<typename V::value_type> temp;
    while (to_visit.size() > 0) {
      const auto [lay, node, child] = to_visit.top();
      to_visit.pop();
      const auto parent = layers[lay][node];
      // first time we see parent? get its label
      if (child == 0)
        temp.push_back(parent.label);

      // base case: reached the bottom layer
      size_t *children = child_buffer + parent.cbuffer_offset;
      if (lay == this->dim - 1) {
        assert(child == 0);
        for (size_t i = 0; i < parent.numchild; i++) {
          auto bottom_node = layers[lay + 1][children[i]];
          temp.push_back(bottom_node.label);
          std::vector<typename V::value_type> cpy{temp};
          res.push_back(V(std::move(cpy)));
          temp.pop_back();
        }
        temp.pop_back(); // done with this parent
      } else {           // recursive case
        // Either we're done with this node and we just mark it as visited or
        // we need to keep it and we add it's next son
        if (child < parent.numchild) {
          to_visit.push({lay, node, child + 1});
          to_visit.push({lay + 1, children[child], 0});
        } else {
          temp.pop_back(); // done with this parent
        }
      }
    }
    return res;
  }

  bool cover_vector(size_t root, V covered) {
    // Stack with tuples (layer, node id, child id)
    std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

    // Add all roots at dimension 0 such that their labels cover the first
    // component of the given vector
    st_node &rootNode = layers[0][root];
    size_t *root_children = child_buffer + rootNode.cbuffer_offset;
    for (size_t i = 0; i < rootNode.numchild; i++) {
      assert(root_children[i] < layers[1].size());
      if (covered[0] <= layers[1][root_children[i]].label)
        to_visit.push({1, root_children[i], 0});
    }

    while (to_visit.size() > 0) {
      const auto [lay, node, child] = to_visit.top();
      to_visit.pop();
      const auto parent = layers[lay][node];
      size_t *children = child_buffer + parent.cbuffer_offset;

      // base case: reached the bottom layer
      if (lay == this->dim - 1) {
        assert(child == 0);
        // it is sufficient to check the last child
        size_t i = parent.numchild - 1;
        auto bottom_node = layers[lay + 1][children[i]];
        if (covered[lay] <= bottom_node.label)
          return true;
      } else { // recursive case
        // Either we're done with this node and we just mark it as visited or
        // we need to keep it and we add it's next son
        size_t c = child;
        if (c == 0) {
          size_t i = parent.numchild - 1;
          auto child_node = layers[lay + 1][children[i]];
          // early exit if the largest child is smaller
          if (covered[lay + 1] > child_node.label)
            continue;
          // otherwise, we find the first index where the order holds
          int left = 0;
          int right = parent.numchild - 2;
          while (left <= right) {
            int mid = left + (right - left) / 2;
            int midVal = layers[lay + 1][children[mid]].label;
            if (covered[lay + 1] <= midVal)
              right = mid - 1;
            else
              left = mid + 1;
          }
          c = left;
        }
        assert(c < parent.numchild);
        to_visit.push({lay, node, c + 1});
        to_visit.push({lay + 1, children[c], 0});
      }
    }
    return false;
  }

  template <std::ranges::input_range R> size_t add_vectors(R &&elements) {
    assert(!layers.empty());

    auto elementVec = std::move(elements);
    // We start a Trie encoded as a map from prefixes to sets of indices of
    // the original vectors, we insert the root too: an empty prefix mapped to
    // the set of all indices
    std::vector<size_t> vectorIds(elementVec.size());
    std::iota(vectorIds.begin(), vectorIds.end(), 0);

    size_t root_id = build_node(vectorIds, 0, elementVec);
    return root_id;
  }
};

template <Vector V>
inline std::ostream &operator<<(std::ostream &os, const sforest<V> &f) {
  for (auto &&el : f.get_all())
    os << el << std::endl;

  return os;
}

} // namespace posets::utils
