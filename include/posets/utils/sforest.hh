#pragma once

#include <cassert>
#include <iostream>
#include <ranges>
#include <stack>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <map>
#include <vector>
#include <optional>

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
    bool isEnd;
    size_t cbuffer_offset;
    size_t numchild;
  };

  struct st_hash {
    sforest* f;
    st_hash (sforest* that) : f{that} {}

    size_t operator()(const st_node &k) const {
      size_t res = std::hash<int>()(k.label);
      size_t* children = f->child_buffer + k.cbuffer_offset;
      for (size_t i = 0; i < k.numchild; i++)
        res ^= std::hash<size_t>()(children[i]) << (i + 1);
      return res;
    }
  };

  struct st_equal {
    sforest* f;
    st_equal (sforest* that) : f{that} {}

    bool operator()(const st_node &lhs, const st_node &rhs) const {
      if (lhs.label != rhs.label or lhs.numchild != rhs.numchild)
        return false;
      size_t* lhs_children = f->child_buffer + lhs.cbuffer_offset;
      size_t* rhs_children = f->child_buffer + rhs.cbuffer_offset;
      for (size_t i = 0; i < lhs.numchild; i++)
        if (lhs_children[i] != rhs_children[i])
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

  std::vector<std::unordered_map<st_node, size_t, st_hash, st_equal>> inverse;

  void init(int k, size_t dim) {
    this->k = k;
    this->dim = dim;
    layers = new st_node *[dim+2];
    layer_size = new size_t[dim+2];
    layer_nxt = new size_t[dim+2];
    for (size_t i = 0; i < dim + 2; i++)
      inverse.emplace_back(INIT_LAYER_SIZE, st_hash(this), st_equal(this));
    for (size_t i = 0; i < dim + 2; i++) {
      layer_size[i] = INIT_LAYER_SIZE;
      layers[i] = new st_node[layer_size[i]];
      layer_nxt[i] = 0;
    }
    cbuffer_size = INIT_LAYER_SIZE * (k + 1);
    child_buffer = new size_t[cbuffer_size];
    cbuffer_nxt = 0;
  }

  

  std::optional<size_t> hasSon(st_node& node, size_t childLayer, int val) {
    if(node.numchild == 0) return std::nullopt;
    
    size_t left = 0;
    size_t right = node.numchild - 1;
    size_t* children = child_buffer + node.cbuffer_offset;

    while (left <= right) {
        size_t mid = left + (right - left) / 2;
        assert (mid < node.numchild);
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

  void addSon(st_node& node, size_t sonLayer, size_t son) {
    // Find the insertion point using binary search
    int left = 0;
    int right = node.numchild - 1;
    st_node& sonNode = layers[sonLayer][son];
    size_t* children = child_buffer + node.cbuffer_offset;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        int midVal = layers[sonLayer][children[mid]].label;
        
        if (sonNode.label == midVal) {
            // This is not supposed to happen -> may add the union approach here directly
            assert (false);
            return;
        }
        else if (sonNode.label < midVal) {
            right = mid - 1;
        }
        else {
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

  size_t addNode(st_node& node, size_t destinationLayer) {
    auto existingNode = inverse[destinationLayer].find(node);
    if(existingNode == inverse[destinationLayer].end()) {
      size_t newID = layer_nxt[destinationLayer];
      inverse[destinationLayer][node] = newID;
      layer_nxt[destinationLayer]++;
      return newID;
    }
    else {
      return existingNode->second;
    }
  }

  size_t addChildren() {
    // Double the child buffer if it is full
    if (cbuffer_nxt + k >= cbuffer_size) {
      size_t* newBuffer = new size_t[cbuffer_size * 2];
      for(size_t i = 0; i < cbuffer_size; i++) {
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

  void removeSon(st_node& node, st_node* layer, size_t son) {
      // Maybe not even necessary? Otherwise need a replaceSon function
  }

  size_t copy(st_node& node, size_t destinationLayer) {
    st_node& newNode = layers[destinationLayer][layer_nxt[destinationLayer]];
    newNode.label = node.label;
    newNode.isEnd = node.isEnd;
    newNode.numchild = node.numchild;
    newNode.cbuffer_offset = addChildren();
    
    size_t* children = child_buffer + node.cbuffer_offset;
    for(size_t i = 0; i < node.numchild; i++) {
      st_node& childNode = layers[destinationLayer + 1][children[i]];
      size_t newSon = copy(childNode, destinationLayer + 1);
      addSon(newNode, destinationLayer + 1, newSon);
    }
    return addNode(newNode, destinationLayer);
  }

  size_t node_union(size_t n_s, size_t n_t, size_t destinationLayer) {
    st_node& node_s = layers[destinationLayer][n_s];
    st_node& node_t = layers[destinationLayer][n_t];
    st_node& newNode = layers[destinationLayer][layer_nxt[destinationLayer]];
    if(node_s.isEnd) {
      newNode.isEnd = true;
    }
    else {
      newNode.label = node_s.label;
      newNode.numchild = node_s.numchild;
      newNode.cbuffer_offset = addChildren();
      size_t s_s{ 0 };
      size_t s_t{ 0 };
      size_t newChild;

      size_t* node_s_children = child_buffer + node_s.cbuffer_offset;
      size_t* node_t_children = child_buffer + node_t.cbuffer_offset;
      while(s_s < node_s.numchild || s_t < node_t.numchild) {
        // Case one: One of the lists is done iterating, copy the nodes
        // without match
        st_node& son_s = layers[destinationLayer + 1][node_s_children[s_s]];
        st_node& son_t = layers[destinationLayer + 1][node_s_children[s_t]];
        if(s_s == node_s.numchild) {
          newChild = copy(son_t, destinationLayer + 1);
          s_t++;
        }
        else if(s_t == node_t.numchild) {
          newChild = copy(son_s, destinationLayer + 1);
          s_s++;
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
        else if (node_s.label > node_t.label) {
          newChild = copy(son_s, destinationLayer + 1);
          s_s++;
        } else {
          newChild = copy(son_t, destinationLayer + 1);
          s_t++;
        }

        addSon(newNode, destinationLayer + 1, newChild);
      }
    }
    return addNode(newNode, destinationLayer);

  }

  std::map<std::vector<size_t>, std::vector<size_t>>
  buildLayer(std::map<std::vector<size_t>, std::vector<size_t>> &layerData,
             int currentLayer, const auto& elementVec) {
    // Init result
    std::map<std::vector<size_t>, std::vector<size_t>> resultNodes{};

    // We have built all layers, so we create the final layer with the EoL
    // node
    if (static_cast<size_t>(currentLayer) == this->dim + 1) {
      assert (layer_nxt[currentLayer] < layer_size[currentLayer]);
      st_node& endNode = layers[currentLayer][layer_nxt[currentLayer]];
      endNode.isEnd = true;
      endNode.label = -1;
      endNode.numchild = 0;
      // FIXME: Don't we need to increase layer_nxt and possible reallocate
      // memory?

      size_t endNodeID = addNode(endNode, currentLayer);
      for (auto const &[n, children] : layerData) {
        resultNodes[n].push_back(endNodeID);
      }

    } else {
      // Re-Partition the nodes based on a prefix with one more entry
      std::map<std::vector<size_t>, std::vector<size_t>> newPartition{};
      for (auto const &[n, children] : layerData) {
        for (auto const &child : children) {
          std::vector<size_t> newPrefix{n};
          newPrefix.push_back(elementVec[child][currentLayer-1]);
          newPartition[newPrefix].push_back(child);
        }
      }

      // Create the next Layer and fill it with the child nodes
      auto childNodes = buildLayer(newPartition, currentLayer+1, elementVec);

      // Create a node for each partition and add the children from that
      // partition
      for (auto const &[n, children] : childNodes) {
        assert (layer_nxt[currentLayer] < layer_size[currentLayer]);
        st_node& newNode = layers[currentLayer][layer_nxt[currentLayer]];
        newNode.label = n.back();
        newNode.numchild = 0;
        newNode.cbuffer_offset = addChildren();

        for (auto const &child : children) {
          // A ST node cannot have two sons with the same value, so we check
          int sonValue = layers[currentLayer +1][child].label;
          std::optional<size_t> sameValueSon = hasSon(newNode, currentLayer + 1, sonValue);
          if (!sameValueSon.has_value()) {
            // Just add the node as son if all is well
            addSon(newNode, currentLayer+1, child);
          } else {
            // If there already is a son with that value, remove it and
            // compute a union-node of the two instead
            size_t unionSon =
                node_union(sameValueSon.value(), child, currentLayer + 1);
            addSon(newNode, currentLayer+1, unionSon);
            // TODO: addSon is probably wrong here, need to replace the original
          }
        }
        // Add the created node to the layer and insert it into the result set
        std::vector<size_t> newPrefix{n};
        newPrefix.pop_back();
        resultNodes[newPrefix].push_back(addNode(newNode, currentLayer));        
      }
    }

    return resultNodes;
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
    delete[] child_buffer;
  }

  std::vector<V> get_all() const {
    // Stack with tuples (layer, node id, child id)
    std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

    // Add all roots at dimension 0
    for (size_t i = 0; i < layer_nxt[1]; i++)
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
      size_t* children = child_buffer + parent.cbuffer_offset;
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
    st_node& rootNode = layers[0][root];
    size_t* root_children = child_buffer + rootNode.cbuffer_offset;
    for (size_t i = 0; i < rootNode.numchild; i++) {
      assert(root_children[i] < layer_nxt[1]);
      if (covered[0] <= layers[1][root_children[i]].label)
        to_visit.push({1, root_children[i], 0});
    }

    while (to_visit.size() > 0) {
      const auto [lay, node, child] = to_visit.top();
      to_visit.pop();
      const auto parent = layers[lay][node];
      size_t* children = child_buffer + parent.cbuffer_offset;

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
          // find the first index where the order holds
          // FIXME: this could be a binary search
          if (covered[lay+1] > child_node.label)
            continue;
          do {
            child_node = layers[lay + 1][children[c]];
            c++;
          } while (covered[lay + 1] > child_node.label);
          c--;
        }
        if (c < parent.numchild) {
          to_visit.push({lay, node, c + 1});
          to_visit.push({lay + 1, children[c], 0});
        }
      }
    }
    return false;
  }

  template <std::ranges::input_range R>
  size_t add_vectors(R&& elements) {
    assert (layers != nullptr);
    st_node* rootLayer = layers[0];
    assert (rootLayer != nullptr);
    st_node& root = rootLayer[layer_nxt[0]];
    root.numchild = 0;
    root.cbuffer_offset = addChildren();

    auto elementVec = std::move(elements);
    // We start a Trie encoded as a map from prefixes to sets of indices of
    // the original vectors, we insert the root too: an empty prefix mapped to
    // the set of all indices
    std::vector<size_t> vectorIds(elementVec.size());
    std::iota(vectorIds.begin(), vectorIds.end(), 0);
    std::vector<size_t> pref{};
    std::map<std::vector<size_t>, std::vector<size_t>> vectorData{
        {pref, vectorIds}};
    // Now, we make a recursive call to build the next layer
    // FIXME: instead of expecting recursive calls to return the children, so
    // that we add them as such here, we can as the recursive call to create
    // the node (if needed) and to adds the children too (this avoid
    // copying/passing around vectors)
    auto children = buildLayer(vectorData, 1, elementVec);
    for (auto const &child : children[pref]) {
      int sonValue = layers[1][child].label;
      std::optional<size_t> sameValueSon = hasSon(root, 1, sonValue);
      if (!sameValueSon.has_value()) {
        // Just add the node as son if all is well
        addSon(root, 1, child);
      } else {
        // If there already is a son with that value, remove it and
        // compute a union-node of the two instead
        size_t unionSon =
            node_union(sameValueSon.value(), child, 1);
        addSon(root, 1, unionSon);
        // TODO: again probably need to replace instead of add
      }
    }
    
    size_t root_id = addNode(root, 0);
    return root_id;
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
