#pragma once

#include <algorithm>
#include <cassert>
#include <map>
#include <numeric>
#include <ranges>
#include <vector>

namespace posets::utils {

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
  size_t *layer_size;
  size_t *layer_nxt;
  size_t *child_buffer;
  size_t *cbuffer_size;
  size_t *cbuffer_nxt;

  std::map<std::vector<size_t>, size_t> cache;

public:
  template <std::ranges::input_range R>
  sharingtree(R &&elements)
      : maxLayers{elements.begin()->size()},
        maxNodes{elements.begin()->size() * elements.size() + 2},
        maxSons{elements.begin()->size() * elements.begin()->size() *
                elements.size()} {
    layerBuffer = new st_layer[maxLayers];
    nodeBuffer = new st_node[maxNodes];
    sonBuffer = new st_son[maxSons];

    // Create root with arbitrary value, will not be read
    root = createNode(0, false, true);
    st_layer_ptr layer1 = addFirstLayer();

    elementVec = std::move(elements);
    std::vector<size_t> vectorIds(this->elementVec.size());
    std::iota(vectorIds.begin(), vectorIds.end(), 0);
    std::vector<size_t> pref{};
    std::map<std::vector<size_t>, std::vector<size_t>> vectorData{
        {pref, vectorIds}};
    auto children = buildLayer(layer1, vectorData, maxLayers);
    for (auto const &child : children[pref]) {
      // A ST node cannot have two sons with the same value, so we check
      st_node_ptr sameValueSon = hasSon(root, child->val);
      if (sameValueSon == nullptr) {
        // Just add the node as son if all is well
        addSon(root, child);
      } else {
        // If there already is a son with that value, remove it and compute a
        // union-node of the two instead
        removeSon(root, sameValueSon);
        st_node_ptr unionSon = node_union(sameValueSon, child, *this, layer1);
        addSon(root, unionSon);
      }
    }
  }

  ~sharingtree() {
    delete[] this->sonBuffer;
    delete[] this->nodeBuffer;
    delete[] this->layerBuffer;
  }

  /*
   * Inclusion
   */
  bool st_includes(V &x) {
    st_son_ptr n = this->root->firstSon;
    while (n != nullptr) {
      // If x is nonempty, but n is already the end-node, we have to check the
      // next branch
      if (x.size() >= 1 && n->node->isEnd) {
        n = n->nextSon;
      } else if (node_includes(n->node, x)) {
        return true;
      } else {
        n = n->nextSon;
      }
    }
    // If we get here, no match was found
    return false;
  }

  /*
   * Union algorithm
   *
   */

  sharingtree st_union(sharingtree &T) {
    sharingtree U{*this, T};
    U.root = U.createNode(0, false, true);

    st_son_ptr n_s = this->root->firstSon;
    st_son_ptr n_t = T.root->firstSon;

    if (n_s != nullptr && n_t != nullptr) {
      st_layer_ptr newLayer = U.addLastLayer();
      st_node_ptr newChild{nullptr};
      while (n_s != nullptr || n_t != nullptr) {
        // Case one: One of the lists is done iterating, copy the nodes
        // without match
        if (n_s == nullptr) {
          newChild = copyIfNotSimulated(n_t->node, U, newLayer, U.root);
          n_t = n_t->nextSon;
        } else if (n_t == nullptr) {
          newChild = copyIfNotSimulated(n_s->node, U, newLayer, U.root);
          n_s = n_s->nextSon;
        }
        // Case two: The values are identical, we union the two nodes
        else if (n_s->node->val == n_t->node->val) {
          newChild = node_union(n_s->node, n_t->node, U, newLayer);
          n_s = n_s->nextSon;
          n_t = n_t->nextSon;
        }
        // Case three: One of the lists is "ahead", we know because the nodes
        // in a layer are ordered
        else if (n_s->node->val > n_t->node->val) {
          newChild = copyIfNotSimulated(n_s->node, U, newLayer, U.root);
          n_s = n_s->nextSon;
        } else {
          newChild = copyIfNotSimulated(n_t->node, U, newLayer, U.root);
          n_t = n_t->nextSon;
        }

        // Add as a  son
        if (newChild != nullptr) {
          U.addSon(U.root, newChild);
        }
      }
    }
    return U;
  }

  /*
   * Intersection algorithm
   *
   */
  sharingtree st_intersect(sharingtree &T) {
    sharingtree I{*this, T};
    I.root = I.createNode(0, false, true);

    st_son_ptr n_s = this->root->firstSon;
    st_layer_ptr newLayer = I.addLastLayer();

    while (n_s != nullptr) {
      st_son_ptr n_t = T.root->firstSon;
      while (n_t != nullptr) {
        st_node_ptr newChild =
            node_intersect(n_s->node, n_t->node, I, newLayer, I.root);
        if (newChild != nullptr) {
          I.addSon(I.root, newChild);
        }
        n_t = n_t->nextSon;
      }
      n_s = n_s->nextSon;
    }

    bool stop = false;
    while (!stop) {
      if (I.lastLayer == nullptr) {
        stop = true;
      } else if (I.lastLayer->firstNode != nullptr) {
        stop = true;
      } else {
        I.deleteLastLayer();
      }
    }

    reduce(I);
    return I;
  }
};
} // namespace posets::utils
