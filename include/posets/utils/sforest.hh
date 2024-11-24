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
  template <Vector V2>
  friend std::ostream &operator<<(std::ostream &os, const utils::sforest<V2> &f);

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

  /*
   Simulation: check if n1 simulates n2
  */
  bool simulates(st_node &n1, st_node &n2, size_t sonLayer) {
    // If the node is in the last layer, we just check the labels
    if(sonLayer == this->dim + 1) {
      return n1.label >= n2.label;
    }
    size_t* n1_children = child_buffer + n1.cbuffer_offset;
    size_t* n2_children = child_buffer + n2.cbuffer_offset;
    if(n1.label < n2.label || 
        layers[sonLayer][n1_children[0]].label < layers[sonLayer][n2_children[n2.numchild - 1]].label)
    {
      // If the label of n1 is too small or its largest child is already smaller than n2's smallest,
      // we already know it can't simulate
      return false;
    }

    // Check if we can find a corresponding son of n1 for every son of n2
    for (size_t s2 = 0; s2 < n2.numchild; s2++)
    {
      bool found = false;
      for (size_t s1 = 0; s1 < n1.numchild; s1++)
      {
        if(simulates(layers[sonLayer][n1_children[s1]], layers[sonLayer][n2_children[s2]], sonLayer + 1)) {
          found = true;
          break;
        }
      }
      // We checked all sons of n1 and there was no match, it can't simulate
      if(!found) return false;
    }

    return true;
  }

  void addSon(st_node &node, size_t sonLayer, size_t son) {
    int last = node.numchild - 1;
    st_node &sonNode = layers[sonLayer][son];
    size_t *children = child_buffer + node.cbuffer_offset;

    // the new node is smaller than all existing ones
    assert(last == -1 or
           layers[sonLayer][children[last]].label > sonNode.label);

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

  size_t addChildren(int numChild) {
    // Double the child buffer if it is full
    if (cbuffer_nxt + numChild >= cbuffer_size) {
      size_t *newBuffer = new size_t[cbuffer_size * 2];
      for (size_t i = 0; i < cbuffer_size; i++) {
        newBuffer[i] = child_buffer[i];
      }
      cbuffer_size *= 2;
      delete[] child_buffer;
      child_buffer = newBuffer;
    }
    size_t res = cbuffer_nxt;
    cbuffer_nxt += numChild;
    return res;
  }

  std::optional<size_t> add_if_not_simulated(st_node &node, size_t destinationLayer, st_node &father) {
    size_t* siblings = child_buffer + father.cbuffer_offset;
    for(size_t s = 0; s < father.numchild; s++) {
      if(simulates(layers[destinationLayer][siblings[s]], node, destinationLayer + 1)) {
        return std::nullopt;
      }
    }

    return addNode(node, destinationLayer);
  }

  size_t node_union(size_t n_s, size_t n_t, size_t destinationLayer) {
    st_node &node_s = layers[destinationLayer][n_s];
    st_node &node_t = layers[destinationLayer][n_t];
    st_node newNode{node_s.label, 0};
    // We haven't reached the bottom of the tree, we need to add children
    if (destinationLayer < this->dim) {
      newNode.cbuffer_offset = addChildren(node_s.numchild + node_t.numchild);
      size_t s_s{0};
      size_t s_t{0};
      size_t newChild;

      size_t *node_s_children = child_buffer + node_s.cbuffer_offset;
      size_t *node_t_children = child_buffer + node_t.cbuffer_offset;
      while (s_s < node_s.numchild || s_t < node_t.numchild) {
        // Case one: One of the lists is done iterating, copy the nodes
        // without match
        st_node &son_s = layers[destinationLayer + 1][node_s_children[s_s]];
        st_node &son_t = layers[destinationLayer + 1][node_t_children[s_t]];
        if (s_s == node_s.numchild) {
          auto copyRes = add_if_not_simulated(son_t, destinationLayer + 1, newNode);
          s_t++;
          if(!copyRes.has_value()) continue;
          newChild = copyRes.value();
        } else if (s_t == node_t.numchild) {
          auto copyRes = add_if_not_simulated(son_s, destinationLayer + 1, newNode);
          s_s++;
          if(!copyRes.has_value()) continue;
          newChild = copyRes.value();
        }
        // Case two: The values are identical, we union the two nodes
        else if (son_s.label == son_t.label) {
          newChild = node_union(node_s_children[s_s], node_t_children[s_t],
                                destinationLayer + 1);
          s_s++;
          s_t++;
        }
        // Case three: One of the lists is "ahead", we know because the nodes
        // in a layer are ordered
        else if (son_s.label > son_t.label) {
          auto copyRes = add_if_not_simulated(son_s, destinationLayer + 1, newNode);
          s_s++;
          if(!copyRes.has_value()) continue;
          newChild = copyRes.value();
        } else {
          auto copyRes = add_if_not_simulated(son_t, destinationLayer + 1, newNode);
          s_t++;
          if(!copyRes.has_value()) continue;
          newChild = copyRes.value();
        }

        addSon(newNode, destinationLayer + 1, newChild);
      }
    }
    return addNode(newNode, destinationLayer);
  }

  std::optional<size_t> node_intersect(size_t n_s, size_t n_t, size_t destinationLayer, std::optional<st_node> father) {
    st_node &node_s = layers[destinationLayer][n_s];
    st_node &node_t = layers[destinationLayer][n_t];
    st_node newNode{std::min(node_s.label, node_t.label), 0};
    // We haven't reached the bottom of the tree, we need to add children
    if(destinationLayer < this->dim) {
      newNode.cbuffer_offset = addChildren(node_s.numchild + node_t.numchild);

      size_t *node_s_children = child_buffer + node_s.cbuffer_offset;
      size_t *node_t_children = child_buffer + node_t.cbuffer_offset;
      for (size_t s_s = 0; s_s < node_s.numchild; s_s++)
      {
        for (size_t s_t = 0; s_t < node_t.numchild; s_t++)
        {
          auto intersectRes = node_intersect(node_s_children[s_s], node_t_children[s_t],
                                destinationLayer + 1, newNode);
          if(intersectRes.has_value()) {
            addSon(newNode, destinationLayer + 1, intersectRes.value());
          }
        }
      }
      
      if(newNode.numchild == 0) {
        return std::nullopt;
      }
    }
    if(father.has_value()) {
      return add_if_not_simulated(newNode, destinationLayer, father.value());
    }
    else {
      return addNode(newNode, destinationLayer);
    }
    
  }

  /* Recursive creation of nodes of Trie while using the inverse map to avoid
   * creating duplicate nodes in terms of (residual/right) language. This
   * results on the creation of the minimal DFA for the set of vectors.
   *
   * Complexity: The implementation below has complexity O(n.d.lg(n)) where n is
   * the number of vectors, d is the number of dimensions. The nd factor is
   * nor surprising: it already corresponds to the set of prefixes of the set
   * of vectors. The lg(n) factor comes from the use of an (ordered) map at
   * every recursive step to (re)partition the vectors for the children based
   * on the list itself. It can be traded off by a factor of k (see TODO in
   * code).
   */
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
      // TODO: Try to do constant-access bucketing based on the value of k.
      // Probably won't pay off unless the set of vectors we are adding is
      // dense in most components.
      std::map<typename V::value_type, std::vector<size_t>, std::greater<typename V::value_type> >
          newPartition{}; // Partition and order the future children
      for (auto const &vec : vecs) {
        newPartition[elementVec[vec][currentLayer]].push_back(vec);
      }
      newNode.cbuffer_offset = addChildren(newPartition.size());

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

  std::vector<V> get_all(std::optional<size_t> root={}) const {
    // Stack with tuples (layer, node id, child id)
    std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

    if (root) {
      to_visit.push({1, root.value(), 0});
    } else {
      // Add all roots at dimension 0
      for (size_t i = 0; i < layers[1].size(); i++)
        to_visit.push({1, i, 0});
    }

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

  void print_children(size_t n, size_t layer) {
    assert(layer <= this->dim);
    st_node& node = layers[layer][n];
    size_t* children = child_buffer + node.cbuffer_offset;
    std::cout << layer << "." << n << " [" << node.label << "] -> ( ";
    for (size_t i = 0; i < node.numchild; i++)
    {
      print_children(children[i], layer + 1);
      if(i < node.numchild -1) {
        std::cout << ", ";
      }
    }
    std::cout << " ) ";
  }

  size_t st_union(size_t root1, size_t root2) {
    return node_union(root1, root2, 0);
  }

  size_t st_intersect(size_t root1, size_t root2) {
    return node_intersect(root1, root2, 0, std::nullopt).value();
  }

  /* Recursive domination check of given vector by vectors in the language of
   * the given tree root.
   *
   * Complexity: The implementation below has complexity O(n.d) where n is
   * the number of vectors, d is the number of dimensions. This is the same as
   * list-based data structures for vectors. However, since it is a DFS of the
   * DFA representation of (bisimulation non-dominated) vectors, it could be
   * exponentially faster than this.
   */
  bool covers_vector(size_t root, V covered) {
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
      assert(lay < this->dim);
      assert(node < layers[lay].size());
      to_visit.pop();
      const auto parent = layers[lay][node];
      assert(child < parent.numchild);
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
        if (c + 1 < parent.numchild)
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

  os << "Layers:" << std::endl;
  for(auto& l : f.layers) {
    for(auto& n: l) {
      os << n.label << " ";
    }
    os << std::endl;
  }
  return os;

}

} // namespace posets::utils
