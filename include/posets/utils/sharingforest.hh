#pragma once

#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <ranges>
#include <stack>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <boost/functional/hash.hpp>

#include <posets/concepts.hh>

namespace posets::utils {

// We will be (micro)managing (C-style) dynamic memory for the set of nodes we
// keep in each layer. The initial number of nodes allowed is determined by
// the number below (multiplied by the number of layers = the dimension + 1).
// When the buffer is full, we double its capacity and copy everything to the
// new block of reserved memory.
size_t INIT_LAYER_SIZE = 100;
unsigned INIT_MAX_CHILDREN = 10;

// Forward definition for the operator<<
template <Vector> class sharingforest;

template <Vector V>
std::ostream &operator<<(std::ostream &os, const utils::sharingforest<V> &f);

template <Vector V> class sharingforest {
private:
  template <Vector V2>
  friend std::ostream &operator<<(std::ostream &os, const utils::sharingforest<V2> &f);

  size_t dim;

  struct st_node {
    typename V::value_type label;
    size_t numchild;
    size_t cbuffer_offset;
  };

  // NOTE: We will be keeping a unique table of st_nodes using an unordered
  // map, for this we need a hash function and an equivalence operation.
  // Since nodes keep offsets instead of pointers, we need both the hash and
  // equivalence operations to be able to access the unique table. Hence, we
  // keep a reference to the table in them.
  struct st_hash {
    sharingforest *f;
    st_hash(sharingforest *that) : f{that} {}

    size_t operator()(const st_node &k) const {
      size_t res = std::hash<typename V::value_type>()(k.label);
      size_t *children = f->child_buffer + k.cbuffer_offset;
      for (size_t i = 0; i < k.numchild; i++)
        res ^= std::hash<size_t>()(children[i]) << (i + 1);
      return res;
    }
  };

  struct st_equal {
    sharingforest *f;
    st_equal(sharingforest *that) : f{that} {}

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

  // This is a cache/hash-map to get in-layer node identifiers from their
  // signature
  std::vector<std::unordered_map<st_node, size_t, st_hash, st_equal>> inverse;
  // Second cache/hash-map to check for a pair of node(-identifiers) in a layer
  // whether there is a simulation relation in the left-to-right direction
  std::vector<std::unordered_map<std::pair<size_t, size_t>, bool, boost::hash<std::pair<size_t, size_t>>>> simulating;

  void init(size_t dim) {
    this->dim = dim;
    layers.resize(dim + 1);

    for (size_t i = 0; i < dim + 1; i++) {
      inverse.emplace_back(INIT_LAYER_SIZE, st_hash(this), st_equal(this));
      simulating.emplace_back();
    }

    cbuffer_size = INIT_LAYER_SIZE * INIT_MAX_CHILDREN;
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
      typename V::value_type midVal = layers[childLayer][children[mid]].label;

      if (midVal == val) {
        return mid;
      } else if (midVal > val) {
        left = mid + 1;
      } else {
        right = mid - 1;
      }
    }

    return std::nullopt;
  }

  bool simulates(size_t n1idx, size_t n2idx, size_t layidx) {
    st_node n1 = layers[layidx][n1idx];
    st_node n2 = layers[layidx][n2idx];
    // If the node is in the last layer, we just check the labels
    if (layidx == this->dim) {
      return n1.label >= n2.label;
    }
    size_t* n1_children = child_buffer + n1.cbuffer_offset;
    size_t* n2_children = child_buffer + n2.cbuffer_offset;
    if (n1.label < n2.label || 
        layers[layidx + 1][n1_children[0]].label <
        layers[layidx + 1][n2_children[n2.numchild - 1]].label) {
      // If the label of n1 is too small or its largest child is already smaller than n2's smallest,
      // we already know it can't simulate
      return false;
    }

    // Stack contains node and child ID of S, node and child ID of T, layer
    std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;
    current_stack.push({n1idx, 0, n2idx, 0, layidx});
 
    while (not current_stack.empty()) {
      auto [n1idx, c1, n2idx, c2, layidx] = current_stack.top();
      current_stack.pop();
      n1 = layers[layidx][n1idx];
      n2 = layers[layidx][n2idx];
#ifndef NDEBUG
      std::cout << "Dimension = " << this->dim << std::endl;
      std::cout << "Comparing " << n1.label << " (" << layidx
                << "." << n1idx << ",c=" << c1 << "/" << n1.numchild
                << ") vs " << n2.label << " ("
                << layidx << "." << n2idx << ",c=" << c2
                << "/" << n2.numchild
                << ")" << std::endl;
#endif

      // This is our base case, no more things to check on the n2 side
      // We now claim success for n2 being simulated by n1 and update the top
      // of the stack so that when we go back up in the tree we have one less
      // branch to check on the n2 side
      if (c2 == n2.numchild  || layidx == this->dim) {
        assert(c2 == n2.numchild);
        simulating[layidx][std::make_pair(n1idx, n2idx)] = true;
        if (not current_stack.empty()) {
          auto [m1idx, d1, m2idx, d2, ell] = current_stack.top();
          current_stack.pop();
          current_stack.push({m1idx, 0, m2idx, d2 + 1, ell});
        }

      // Another base case: we've iterated through all the children on the
      // n1 side and failed to find a simulating one
      } else if (c1 == n1.numchild) {
        simulating[layidx][std::make_pair(n1idx, n2idx)] = false;
        return false;

      // The last case is that we have two valid children indices,
      // then we have to go deeper in the product tree (lest the cache saves
      // us, or we find that the subtree will just certainly not satisfy
      // simulation)
      } else {
        n1_children = child_buffer + n1.cbuffer_offset;
        n2_children = child_buffer + n2.cbuffer_offset;
        auto node_pair = std::make_pair(n1_children[c1],
                                        n2_children[c2]);
        auto cached = simulating[layidx + 1].find(node_pair);
        // Did we get lucky with the cache? then push back an updated node
        // with less obligations or keep searching on the n1 side
        if (cached != simulating[layidx + 1].end()) {
          if (cached->second)
            current_stack.push({n1idx, 0, n2idx, c2 + 1, layidx});
          else 
            current_stack.push({n1idx, c1 + 1, n2idx, c2, layidx});
        } else {
          st_node c1_node = layers[layidx + 1][n1_children[c1]];
          st_node c2_node = layers[layidx + 1][n2_children[c2]];
          size_t *c1_children = child_buffer + c1_node.cbuffer_offset;
          size_t *c2_children = child_buffer + c2_node.cbuffer_offset;
          // In some cases, we can eliminate the subtree altogether
          if (c1_node.label < c2_node.label or
              (layidx + 1 < this->dim and
               layers[layidx + 2][c1_children[0]].label <
               layers[layidx + 2][c2_children[c2_node.numchild - 1]].label)) {
            current_stack.push({n1idx, c1 + 1, n2idx, c2, layidx});
          // Alright, we have to go deeper now; and push back what we just
          // popped
          } else {
            current_stack.push({n1idx, c1, n2idx, c2, layidx});
            current_stack.push({n1_children[c1], 0,
                                n2_children[c2], 0,
                                layidx + 1});
          }
        }
      }
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

  void addSonUnordered(st_node &node, size_t sonLayer, size_t son) {
    // Find the insertion point using binary search
    int left = 0;
    int right = node.numchild - 1;
    st_node &sonNode = layers[sonLayer][son];
    size_t *children = child_buffer + node.cbuffer_offset;
    while (left <= right) {
      int mid = left + (right - left) / 2;
      int midVal = layers[sonLayer][children[mid]].label;

      if (sonNode.label == midVal) {
        size_t newSon = node_union(son, children[mid], sonLayer);
        children[mid] = newSon;
        return;
      } else if (midVal < sonNode.label) {
        right = mid - 1;
      } else {
        left = mid + 1;
      }

    }
    // Shift elements in the child buffer to make room for the new child
    for (int i = node.numchild; i > left; i--) {
      children[i] = children[i - 1];
    }

    // Insert the new value
    children[left] = son;
    node.numchild++;
  }

  void addSonIfNotSimulated(size_t node, size_t sonLayer, st_node &father) {
    size_t* siblings = child_buffer + father.cbuffer_offset;
    for(size_t s = 0; s < father.numchild; s++) {
      if (simulates(siblings[s], node, sonLayer)) {
        return;
      }
    }
    addSon(father, sonLayer, node);
    return;
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
      std::memcpy(newBuffer, child_buffer, cbuffer_size * sizeof(size_t));
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
    size_t res = addNode(node, destinationLayer);
    for(size_t s = 0; s < father.numchild; s++) {
      if (simulates(siblings[s], res, destinationLayer)) {
        return std::nullopt;
      }
    }
    return res;
  }

  size_t node_union(size_t ns, size_t nt, size_t destinationLayer) {
    // Stack contains node and child ID of S, node and child ID of T, layer
    std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;

    st_node rootNode1 = layers[destinationLayer][ns];
    st_node rootNode2 = layers[destinationLayer][nt];
    assert(rootNode1.numchild > 0 and rootNode2.numchild > 0);

    current_stack.push({ns, 0, nt, 0, destinationLayer});

    while(not current_stack.empty())
    {
      auto [n_s, c_s, n_t, c_t, layer] = current_stack.top();
      current_stack.pop();
      assert(n_s < layers[layer].size());
      st_node node_s = layers[layer][n_s];
      assert(n_t < layers[layer].size());
      st_node node_t = layers[layer][n_t];

      // It's the first time we see this combination of nodes, so let's insert a draft
      // of it, still dirty, will be cleaned later
      if (c_s == 0 and c_t == 0) {
        assert (node_s.label == node_t.label);
        if (layer < this->dim) {
          layers[layer].emplace_back(node_s.label, 0,
                                     addChildren(node_s.numchild + node_t.numchild));
        } else {
          layers[layer].emplace_back(node_s.label, 0);
        }
      }

      auto &under_construction = layers[layer].back();
      // Base case: We are done with this combo and can clean up
      if (c_s == node_s.numchild and c_t == node_t.numchild) {
        if (layer == destinationLayer) {
          // We skip this step, we clean the root after the loop, but we are done
          break;
        }
        // The very first thing to do is to check whether our draft of node is
        // not a repetition of something in the table already! 
        auto &father = layers[layer - 1].back();
        layers[layer].pop_back();
        auto unionRes = add_if_not_simulated(under_construction,
                                                 layer, father);
        if (unionRes.has_value()) {
          addSon(father, layer, unionRes.value());
        }
      // Recursive step: Either just add the son to the draft node and 
      // continue with the next node, or continue down the tree if 
      // necessary
      } else if (layer < this->dim) { 
        size_t *node_s_children = child_buffer + node_s.cbuffer_offset;
        size_t *node_t_children = child_buffer + node_t.cbuffer_offset;

        // Case 1: One node is "done"
        // We add the existing node (including all its sons!) to the node
        // currently being constructed
        if (c_s == node_s.numchild) {
          addSonIfNotSimulated(node_t_children[c_t], layer + 1, under_construction);
          if (c_t < node_t.numchild) current_stack.push({n_s, c_s, n_t, c_t + 1, layer});
        }
        else if (c_t == node_t.numchild) {
          addSonIfNotSimulated(node_s_children[c_s], layer + 1, under_construction);
          if (c_s < node_s.numchild) current_stack.push({n_s, c_s + 1, n_t, c_t, layer});
        }
        else 
        {
          st_node &son_s = layers[layer + 1][node_s_children[c_s]];
          st_node &son_t = layers[layer + 1][node_t_children[c_t]];
          // Case 2: One child is "ahead", We add just as when one list is done
          if (son_t.label > son_s.label) {
            addSonIfNotSimulated(node_t_children[c_t], layer + 1, under_construction);
            if (c_t < node_t.numchild) current_stack.push({n_s, c_s, n_t, c_t + 1, layer});
          }
          else if (son_s.label > son_t.label) {
            addSonIfNotSimulated(node_s_children[c_s], layer + 1, under_construction);
            if (c_s < node_s.numchild) current_stack.push({n_s, c_s + 1, n_t, c_t, layer});
          }
          // Case 3: The values of the sons match, so we need to continue the recursion
          else {
            assert (son_s.label == son_t.label);
            if (c_s < node_s.numchild and c_t < node_t.numchild) current_stack.push({n_s, c_s + 1, n_t, c_t + 1, layer});
            current_stack.push({node_s_children[c_s], 0, node_t_children[c_t], 0, layer + 1});
          }
        }
      }
    }
    
    // Clean up the root, we remove it for and re-add it, so it is checked whether an identical node exists
    auto constructed_node = layers[destinationLayer].back();
    layers[destinationLayer].pop_back();
    return addNode(constructed_node, destinationLayer);
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
            // Can happen that we insert the same value twice, so we need to check
            auto existingSon = hasSon(newNode, destinationLayer + 1, layers[destinationLayer + 1][intersectRes.value()].label);
            if(existingSon.has_value()) {
              size_t newSon = node_union(existingSon.value(), intersectRes.value(), destinationLayer + 1);
              size_t* newNode_children = child_buffer + newNode.cbuffer_offset;
              newNode_children[existingSon.value()] = newSon;
            }
            else {
              addSon(newNode, destinationLayer + 1, intersectRes.value());
            }
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
    typename V::value_type label{currentLayer == 0
                  ? static_cast<typename V::value_type>(-1)
                  : elementVec[vecs[0]][currentLayer - 1]};
    st_node newNode{label, 0};
    // We have not reached the last layer - so add children
    if (currentLayer < this->dim) {
      if (vecs.size() == 1) {
        auto& vec = elementVec[vecs[0]];
        st_node lastSonNode{vec[vec.size() - 1], 0};
        size_t nextSon = addNode(lastSonNode, this->dim);
        for (size_t i = vec.size() - 1; i > currentLayer; i--) {
          st_node sonNode{vec[i - 1], 0};
          sonNode.cbuffer_offset = addChildren(1);
          addSon(sonNode, i + 1, nextSon);
          nextSon = addNode(sonNode, i);
        }
        newNode.cbuffer_offset = addChildren(1);
        addSon(newNode, currentLayer + 1, nextSon);
      } else {
        // Partition and order the future children
        // TODO: Try to do constant-access bucketing based on the value of k.
        // Probably won't pay off unless the set of vectors we are adding is
        // dense in most components.
        std::map<typename V::value_type,
                 std::vector<size_t>,
                 std::greater<typename V::value_type>> newPartition{};
        for (auto const &vec : vecs) {
          newPartition[elementVec[vec][currentLayer]].push_back(vec);
        }
        newNode.cbuffer_offset = addChildren(newPartition.size());

        for (auto &[n, children] : newPartition) {
          // Build a new son for each individual value at currentLayer + 1
          size_t newSon = build_node(children, currentLayer + 1, elementVec);
          bool found = false;
          size_t* currentChildren = child_buffer + newNode.cbuffer_offset;
          for (size_t s = 0; s < newNode.numchild; s++) {
            if (simulates(currentChildren[s], newSon, currentLayer + 1)) {
              found = true;
              break;
            }
          }
          if (!found)
            addSon(newNode, currentLayer + 1, newSon);
        }
      }
    }
    return addNode(newNode, currentLayer);
  }

public:
  sharingforest() {}

  sharingforest(size_t dim) {
    assert(dim >= 2);
    this->init(dim);
  }

  ~sharingforest() {
    if (layers.empty())
      return;

    delete[] child_buffer;
  }

  std::vector<V> get_all(std::optional<size_t> root={}) const {
    // Stack with tuples (layer, node id, child id)
    std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

    if (root) {
      const st_node& rootNode = layers[0][root.value()];
      size_t *firstChildren = child_buffer + rootNode.cbuffer_offset;
      for (size_t c = 0; c < rootNode.numchild; c++)
      {
        to_visit.push({1, firstChildren[c], 0});
      }
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
#ifndef NDEBUG
    std::cout << std::string(layer, '\t') << layer << "." << n << " [" << static_cast<int>(node.label) << "] -> (" << (layer == this->dim? "" : "\n");
#endif
    for (size_t i = 0; i < node.numchild; i++)
    {
      print_children(children[i], layer + 1);
    }
#ifndef NDEBUG
    std::cout << (layer == this->dim? " " : std::string(layer, '\t')) << " )\n";
#endif
  }

  size_t st_union(size_t root1, size_t root2) {
    return node_union(root1, root2, 0);
  }

  size_t st_intersect(size_t root1, size_t root2) {
    // Stack contains node and child ID of S, node and child ID of T, layer
    std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;

    st_node rootNode1 = layers[0][root1];
    st_node rootNode2 = layers[0][root2];
    assert(rootNode1.numchild > 0 and rootNode2.numchild > 0);
    size_t *root1_children = child_buffer + rootNode1.cbuffer_offset;
    size_t *root2_children = child_buffer + rootNode2.cbuffer_offset;

    // Insert a draft of root node, dirty one without checking if it's there,
    // we'll clean later
    layers[0].emplace_back(static_cast<typename V::value_type>(-1), 0,
                           addChildren(rootNode1.numchild + rootNode2.numchild));

    // We are ready to start a stack-simulated DFS of the synchronized-product
    // of the trees
    for (size_t c_1 = 1; c_1 <= rootNode1.numchild; c_1++) {
      for (size_t c_2 = 1; c_2 <= rootNode2.numchild; c_2++) {
        current_stack.push({root1_children[rootNode1.numchild - c_1], 0,
                            root2_children[rootNode2.numchild - c_2], 0, 1});
      }
    }
 
    while (not current_stack.empty()) {
      auto [n_s, c_s, n_t, c_t, layer] = current_stack.top();
      current_stack.pop();
      assert(n_s < layers[layer].size());
      st_node node_s = layers[layer][n_s];
      assert(n_t < layers[layer].size());
      st_node node_t = layers[layer][n_t];

      // It's the first time we see this product node, so let's insert a draft
      // of it, again dirty
      if (c_s == 0 and c_t == 0) {
        if (layer < this->dim) {
          layers[layer].emplace_back(std::min(node_s.label, node_t.label), 0,
                                     addChildren(node_s.numchild + node_t.numchild));
        } else {
          layers[layer].emplace_back(std::min(node_s.label, node_t.label), 0);
        }
      }

      // This is our base case, we've iterated through all the children in the
      // product tree, it's finally time to clean up stuff
      if (c_s == node_s.numchild and c_t == node_t.numchild) {
        // The very first thing to do is to check whether our draft of node is
        // not a repetition of something in the table already!
        auto under_construction = layers[layer].back();
        auto &father = layers[layer - 1].back();
        layers[layer].pop_back();
        auto intersectRes = add_if_not_simulated(under_construction,
                                                 layer, father);
        if (intersectRes.has_value()) {
          addSonUnordered(father, layer, intersectRes.value());
        }
      // Below we have the "recursive" step in which we put two elements into
      // the stack to remember what was the next child of the node at this
      // level and to go deeper in the tree if needed
      } else if (layer < this->dim) {
        size_t *node_s_children = child_buffer + node_s.cbuffer_offset;
        size_t *node_t_children = child_buffer + node_t.cbuffer_offset;

        if (c_t < node_t.numchild) {
          current_stack.push({n_s, c_s, n_t, c_t + 1, layer});
        } else if (c_s < node_s.numchild) {
          current_stack.push({n_s, c_s + 1, n_t, 0, layer});
        }
        if (c_t < node_t.numchild and c_s < node_s.numchild) {
          current_stack.push({node_s_children[c_s], 0, node_t_children[c_t], 0, layer + 1});
        }
      }
    }

    // Clean up the root too! It's added to the layer, so we just need to make
    // sure it's not there already and remove this second copy if it is
    auto under_construction = layers[0].back();
    layers[0].pop_back();
    return addNode(under_construction, 0);
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
  bool covers_vector(size_t root, const V &covered) {
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
        // it is sufficient to check the first child
        auto bottom_node = layers[lay + 1][children[0]];
        if (covered[lay] <= bottom_node.label)
          return true;
      } else { // recursive case
        // Either we're done with this node and we just mark it as visited or
        // we need to keep it and we add it's next son
        size_t c = child;
        auto child_node = layers[lay + 1][children[c]];
        // early exit if the largest child is smaller
        if (covered[lay] > child_node.label)
          continue;
        
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
#ifndef NDEBUG
    size_t maxlayer = 0;
    size_t totlayer = 0;
    for (size_t i = 0; i < this->dim + 1; i++) {
      totlayer += this->layers[i].size();
      if (maxlayer < this->layers[i].size())
        maxlayer = this->layers[i].size();
    }
    std::cout << "[Forest stats] Max layer size=" << maxlayer
              << " total layers' size=" << totlayer
              << " (bytes per el=" << sizeof(st_node) << ")"
              << std::endl;
    std::cout << "[Forest stats] Child buff size=" << this->cbuffer_size
              << " (bytes per el=" << sizeof(size_t) << ")"
              << std::endl;
#endif
    return root_id;
  }

  /*
    For testing: Check that all children are ordered descending
  */
  bool check_child_order() {
    size_t layer_num = 0;
    for(auto& l : layers) {
      for(st_node& n: l) {
        size_t* children = child_buffer + n.cbuffer_offset;
        for (size_t i = 1; i < n.numchild; i++)
        {
          // There is a child with a larger label than the previous child
          if(layers[layer_num + 1][children[i]].label > layers[layer_num + 1][children[i - 1]].label) {
            return false;
          }
        }
      }
      layer_num++;
    }
    return true;
  }

  /*
    For testing: Check that one root simulates another
  */
  bool check_simulation(size_t n1, size_t n2) {
    return simulates(n1, n2, 0);
  }
};

template <Vector V>
inline std::ostream &operator<<(std::ostream &os, const sharingforest<V> &f) {
  for (auto &&el : f.get_all())
    os << el << std::endl;

  os << "Layers:" << std::endl;
  for(auto& l : f.layers) {
    for(auto& n: l) {
      os << static_cast<int>(n.label) << " ";
    }
    os << std::endl;
  }
  return os;

}

} // namespace posets::utils
