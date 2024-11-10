#pragma once

#include <algorithm>
#include <cassert>
#include <map>
#include <numeric>
#include <ranges>
#include <vector>

#include <posets/concepts.hh>

namespace posets::utils {

template <Vector V> class sharingtree {
private:
  struct st_layer;
  using st_layer_ptr = st_layer *;
  struct st_node;
  using st_node_ptr = st_node *;
  struct st_son;
  using st_son_ptr = st_son *;

  // Variables for management of components
  size_t maxLayers;
  st_layer_ptr layerBuffer;
  int layerBufferIndex{0};
  size_t maxNodes;
  st_node_ptr nodeBuffer;
  int nodeBufferIndex{0};
  size_t maxSons;
  st_son_ptr sonBuffer;
  int sonBufferIndex{0};

  std::vector<V> elementVec;

  st_layer_ptr firstLayer{nullptr};
  st_layer_ptr lastLayer{nullptr};
  st_node_ptr root{nullptr};

  struct st_layer {
    st_node_ptr firstNode{nullptr};
    st_node_ptr lastNode{nullptr};

    st_layer_ptr previousLayer{nullptr};
    st_layer_ptr nextLayer{nullptr};
  };

  struct st_node {
    st_son_ptr firstSon{nullptr};
    size_t val;
    bool isRoot{false};
    bool isEnd{false};

    st_node_ptr nextNode{nullptr};
    // still need the auxiliary values
  };

  struct st_son {
    st_node_ptr node{nullptr};
    st_son_ptr nextSon{nullptr};
  };

  st_node_ptr createNode(size_t val, bool isEnd = false, bool isRoot = false) {
    st_node_ptr newNode = &(nodeBuffer[nodeBufferIndex]);
    nodeBufferIndex++;
    newNode->val = val;
    newNode->isEnd = isEnd;
    newNode->isRoot = isRoot;

    return newNode;
  }

  /*
  Adds the specified node as child
  */
  void addSon(st_node_ptr node, st_node_ptr son) {
    assert(node != nullptr);
    st_son_ptr currentSon = node->firstSon;
    // Check if the new node needs to be the first son
    if (currentSon == nullptr || son->val > currentSon->node->val) {
      st_son_ptr newSon = &(sonBuffer[sonBufferIndex]);
      sonBufferIndex++;
      newSon->node = son;
      newSon->nextSon = currentSon;
      node->firstSon = newSon;
    } else {
      while (currentSon->nextSon != nullptr &&
             currentSon->node->val > son->val) {
        currentSon = currentSon->nextSon;
      }
      // Constraint on successors may not satisfied
      assert(currentSon->node->val != son->val);
      // FIXME: Throw error or insert children as
      // children of node?

      st_son_ptr newSon = &(sonBuffer[sonBufferIndex]);
      sonBufferIndex++;
      newSon->node = son;
      newSon->nextSon = currentSon->nextSon;
      currentSon->nextSon = newSon;
    }
  }

  void removeSon(st_node_ptr node, st_node_ptr son) {
    assert(node != nullptr);
    st_son_ptr currentSon = node->firstSon;
    // Special case: node to remove is the beginning of the successors
    if (son == currentSon->node) {
      node->firstSon = currentSon->nextSon;
    } else {
      while (currentSon->nextSon != nullptr &&
             currentSon->nextSon->node->val > son->val) {
        currentSon = currentSon->nextSon;
      }
      if (currentSon->nextSon->node == son) {
        currentSon->nextSon = currentSon->nextSon->nextSon;
      }
    }
  }

  // FIXME: This implements a linear search on an ORDERED list of integers
  // the children are not necessarily adjacent in memory so a binary
  // search seems (a priori) hard to implement + we know the number of
  // children will be small = so this should not be too bad
  st_node_ptr hasSon(st_node_ptr node, size_t val) {
    assert(node != nullptr);
    st_son_ptr currentSon = node->firstSon;
    while (currentSon != nullptr && currentSon->node->val > val) {
      currentSon = currentSon->nextSon;
    }
    if (currentSon != nullptr && currentSon->node->val == val) {
      return currentSon->node;
    } else {
      return nullptr;
    }
  }

  bool sameSons(st_node_ptr n1, st_node_ptr n2) {
    assert(n1 != nullptr);
    assert(n2 != nullptr);
    st_son_ptr currentSon1 = n1->firstSon;
    st_son_ptr currentSon2 = n2->firstSon;
    if (currentSon1 == nullptr || currentSon2 == nullptr) {
      return currentSon1 == currentSon2;
    }
    bool res = false;
    // Iterate while the node values match
    while (currentSon1->node == currentSon2->node) {
      if (currentSon1->nextSon == nullptr || currentSon2->nextSon == nullptr) {
        // If one of the lists has reached its end, we check whether they are
        // equal (should both be nullptr then) Res only becomes true if both
        // nodes have reached their last son!
        res = currentSon1->nextSon == currentSon2->nextSon;
        break;
      };
      currentSon1 = currentSon1->nextSon;
      currentSon2 = currentSon2->nextSon;
    }
    return res;
  }

  void removeNode(st_layer_ptr layer, st_node_ptr node) {
    assert(layer != nullptr);
    st_node_ptr currentNode = layer->firstNode;
    // Special case: node to remove is the beginning of the layer
    if (currentNode != nullptr && node == currentNode) {
      layer->firstNode = currentNode->nextNode;
    } else if (currentNode != nullptr) {
      while (currentNode->nextNode != nullptr &&
             currentNode->nextNode->val > node->val) {
        currentNode = currentNode->nextNode;
      }
      if (currentNode->nextNode == node) {
        if (currentNode->nextNode == layer->lastNode) {
          layer->lastNode = currentNode;
        } else {
          currentNode->nextNode = currentNode->nextNode->nextNode;
        }
      }
    }
  }

  st_node_ptr addNode(st_layer_ptr layer, st_node_ptr node) {
    assert(layer != nullptr);
    st_node_ptr currentNode = layer->firstNode;
    // check if new node needs to be the first in the layer
    if (currentNode == nullptr || node->val > currentNode->val) {
      node->nextNode = currentNode;
      layer->firstNode = node;
      return node;
    }

    while (currentNode->nextNode != nullptr && currentNode->val > node->val) {
      currentNode = currentNode->nextNode;
    }
    // Check if node with same value and successors as node exists in layer
    while (currentNode->val == node->val) {
      if (sameSons(currentNode, node)) {
        // if it does, the result of the insertion is the existing node
        return currentNode;
      } else {
        // else we have to check if there is another node with the same value
        // and do the same check
        currentNode = currentNode->nextNode;
      }
    }
    // if there exists no matching node, we insert
    node->nextNode = currentNode->nextNode;
    currentNode->nextNode = node;
    // Change last node if the node previously was the last
    if (currentNode == layer->lastNode) {
      layer->lastNode = node;
    }
    return node;
  }

  st_layer_ptr addFirstLayer() {
    st_layer_ptr newFirst = &(layerBuffer[layerBufferIndex]);
    layerBufferIndex++;
    newFirst->nextLayer = firstLayer;
    firstLayer = newFirst;
    return newFirst;
  }

  st_layer_ptr addLastLayer() {
    st_layer_ptr newLast = &(layerBuffer[layerBufferIndex]);
    layerBufferIndex++;
    newLast->previousLayer = lastLayer;
    lastLayer = newLast;
    return newLast;
  }

  void deleteLastLayer() {
    st_layer_ptr newLast = lastLayer->previousLayer;
    newLast->nextLayer = nullptr;
    // delete lastLayer; TODO: need to replace with new mechanism?
    lastLayer = newLast;
  }

  // This is a simple standard copy without checks, we might only use the
  // simulation based copy and this method will be obsolete!
  st_node_ptr copy(st_node_ptr node, sharingtree *S,
                   st_layer_ptr destinationLayer) {
    st_node_ptr newNode = S.createNode(node->val);
    if (node->nextNode != nullptr) {
      st_layer_ptr nextLayer = destinationLayer->nextLayer;
      if (nextLayer == nullptr) {
        nextLayer = S->addLastLayer();
      }
      st_son_ptr son = node->firstSon;
      while (son != nullptr) {
        S.addSon(newNode, copy(son->node, S, nextLayer));
        son = son->nextSon;
      }
    }
    return addNode(destinationLayer, newNode);
  }

  st_node_ptr copyIfNotSimulated(st_node_ptr node, sharingtree &S,
                                 st_layer_ptr destinationLayer,
                                 st_node_ptr father) {
    // Check if there is a node that simulates the node we want, if so, just
    // return that node instead
    assert(father != nullptr);
    st_son_ptr checkNode = father->firstSon;
    while (checkNode != nullptr && checkNode->node->val > node->val) {
      if (simulates(checkNode->node, node)) {
        return checkNode->node;
      } else {
        checkNode = checkNode->nextSon;
      }
    }

    // No simulating node was found, we create a new node and copy its
    // children (if they are not simulated!)
    st_node_ptr newNode = S.createNode(node->val);
    if (node->nextNode != nullptr) {
      st_layer_ptr nextLayer = destinationLayer->nextLayer;
      if (nextLayer == nullptr) {
        nextLayer = S.addLastLayer();
      }
      st_son_ptr son = node->firstSon;
      while (son != nullptr) {
        S.addSon(newNode, copyIfNotSimulated(son->node, S, nextLayer, newNode));
        son = son->nextSon;
      }
    }
    return addNode(destinationLayer, newNode);
  }

  /*
   * Simulation relation check
   */

  bool simulates(st_node_ptr n1, st_node_ptr n2) {
    if (n1->isEnd)
      return n2->isEnd;
    // Preliminary check: Condition 1, value must be greater or equal
    if (n1->val < n2->val)
      return false;
    st_son_ptr s1 = n1->firstSon;
    st_son_ptr s2 =
        n2->firstSon; // Maybe add a check that they aren't just both empty?
                      // Can't be empty bc of ST definition!
    // Check every son of n1
    while (s2 != nullptr) {
      // There was no son of n1 that simulates the son of n2, so we know the
      // node doesn't simulate
      if (s1 == nullptr)
        return false;

      if (simulates(s1->node, s2->node)) {
        // If it simulates, we move on to the next child
        s2 = s2->nextSon;
        // We return to the first child (might not be necessary)
        s1 = n1->firstSon;
      } else {
        // Try the next
        s1 = s1->nextSon;
      }
    }
    // There was a simulating node found for every son, so the node simulates
    return true;
  }

  void node_reduce(st_node_ptr n, st_node_ptr father) {
    if (!n->isEnd) {
      st_son_ptr checkNode = father->firstSon;
      while (checkNode != nullptr && checkNode->node->val > n->val) {
        if (simulates(checkNode->node, n)) {
          removeSon(father, n);
          return;
        } else {
          checkNode = checkNode->nextSon;
        }
      }

      // We did not find a simulating sibling, so we move down one layer
      st_son_ptr s = n->firstSon;
      while (s != nullptr) {
        node_reduce(s->node, n);
      }
    }
  }

  void reduce(sharingtree &S) {
    st_node_ptr n = S.root->firstSon->node;
    while (n != nullptr) {
      node_reduce(n, S.root);
    }
  }

  bool node_includes(st_node_ptr n, const V &x) {
    if (n->val >= *x.begin()) {
      // The current node is a candidate, we check the children
      st_son_ptr s = n->firstSon;
      while (s != nullptr) {
        if (x.size() == 1 && s->node->isEnd) {
          // We are at the end of the input vector and reached an EoL node
          return true;
        } else if (node_includes(s->node, *x.end())) {
          return true;
        } else {
          s = s->nextSon;
        }
      }
    }
    // If we get here, the current branch can be discarded
    return false;
  }

  st_node_ptr node_intersect(st_node_ptr n_s, st_node_ptr n_t, sharingtree &I,
                             st_layer_ptr newLayer, st_node_ptr father) {
    st_node_ptr newNode;
    // Start by inserting EoL if the value is EoL
    if (n_s->isEnd || n_t->isEnd) {
      newNode = I.createNode(0, true);
      I.addNode(newLayer, newNode);
    } else {
      newNode = I.createNode(std::min(n_s->val, n_t->val));
      st_layer_ptr nextLayer = newLayer->nextLayer;
      if (nextLayer == nullptr) {
        nextLayer = I.addLastLayer();
      }
      st_son_ptr s_s = n_s->firstSon;

      while (s_s != nullptr) {
        st_son_ptr s_t = n_t->firstSon;
        while (s_t != nullptr) {
          st_node_ptr newChild =
              node_intersect(s_s->node, s_t->node, I, nextLayer, newNode);
          if (newChild != nullptr) {
            I.addSon(newNode, newChild);
          }
          s_t = s_t->nextSon;
        }
        s_s = s_s->nextSon;
      }

      if (newNode->firstSon != nullptr) {
        // Add if not simulated
        st_son_ptr checkNode = father->firstSon;
        while (checkNode != nullptr && checkNode->node->val > newNode->val) {
          if (simulates(checkNode->node, newNode)) {
            // Discard the node we built and return the simulating one instead
            return checkNode->node;
          } else {
            checkNode = checkNode->nextSon;
          }
        }
        addNode(newLayer, newNode);
      } else {
        newNode = nullptr;
      }
    }
    return newNode;
  }

  st_node_ptr node_union(st_node_ptr n_s, st_node_ptr n_t, sharingtree &U,
                         st_layer_ptr newLayer) {
    st_node_ptr newNode{};
    // Start by inserting EoL if the value is EoL
    if (n_s == nullptr) {
      // TODO
    } else {
      newNode = createNode(n_s->val);
      st_layer_ptr nextLayer = newLayer->nextLayer;
      if (nextLayer == nullptr) {
        nextLayer = U.addLastLayer();
      }
      st_son_ptr s_s = n_s->firstSon;
      st_son_ptr s_t = n_t->firstSon;

      st_node_ptr newChild{nullptr};
      while (s_s != nullptr || s_t != nullptr) {
        // Case one: One of the lists is done iterating, copy the nodes
        // without match
        if (s_s == nullptr) {
          newChild = copyIfNotSimulated(s_t->node, U, nextLayer, newNode);
          s_t = s_t->nextSon;
        } else if (s_t == nullptr) {
          newChild = copyIfNotSimulated(s_s->node, U, nextLayer, newNode);
          s_s = s_s->nextSon;
        }
        // Case two: The values are identical, we union the two nodes
        else if (s_s->node->val == s_t->node->val) {
          newChild = node_union(s_s->node, s_t->node, U, nextLayer);
          s_s = s_s->nextSon;
          s_t = s_t->nextSon;
        }
        // Case three: One of the lists is "ahead", we know because the nodes
        // in a layer are ordered
        else if (s_s->node->val > s_t->node->val) {
          newChild = copyIfNotSimulated(s_s->node, U, nextLayer, newNode);
          s_s = s_s->nextSon;
        } else {
          newChild = copyIfNotSimulated(s_t->node, U, nextLayer, newNode);
          s_t = s_t->nextSon;
        }

        // Add as a  son
        if (newChild != nullptr) {
          U.addSon(newNode, newChild);
        }
      }
      addNode(newLayer, newNode);
    }
    return newNode;
  }

  std::map<std::vector<size_t>, std::vector<st_node_ptr>>
  buildLayer(st_layer_ptr layerK,
             std::map<std::vector<size_t>, std::vector<size_t>> &layerData,
             int k) {
    // Init result
    std::map<std::vector<size_t>, std::vector<st_node_ptr>> resultNodes{};

    // We have built all layers, so we create the final layer with the EoL
    // node
    if (k == 0) {
      st_node_ptr endNode = createNode(0, true);
      endNode = addNode(layerK, endNode);

      for (auto const &[n, children] : layerData) {
        resultNodes[n].push_back(endNode);
      }

    } else {
      // Re-Partition the nodes based on a prefix with one more entry
      std::map<std::vector<size_t>, std::vector<size_t>> newPartition{};
      for (auto const &[n, children] : layerData) {
        for (auto const &child : children) {
          std::vector<size_t> newPrefix{n};
          newPrefix.push_back(elementVec[child][maxLayers - k]);
          newPartition[newPrefix].push_back(child);
        }
      }

      // Create the next Layer and fill it with the child nodes
      st_layer_ptr nextLayer = addLastLayer();
      auto childNodes = buildLayer(nextLayer, newPartition, k - 1);

      // Create a node for each partition and add the children from that
      // partition
      for (auto const &[n, children] : childNodes) {
        st_node_ptr newNode = createNode(n.back());

        for (auto const &child : children) {
          // A ST node cannot have two sons with the same value, so we check
          st_node_ptr sameValueSon = hasSon(newNode, child->val);
          if (sameValueSon == nullptr) {
            // Just add the node as son if all is well
            addSon(newNode, child);
          } else {
            // If there already is a son with that value, remove it and
            // compute a union-node of the two instead
            removeSon(newNode, sameValueSon);
            st_node_ptr unionSon =
                node_union(sameValueSon, child, *this, nextLayer);
            addSon(newNode, unionSon);
          }
        }
        // Add the created node to the layer and insert it into the result set
        newNode = addNode(layerK, newNode);
        std::vector<size_t> newPrefix{n};
        newPrefix.pop_back();
        resultNodes[newPrefix].push_back(newNode);
      }
    }

    return resultNodes;
  }

public:
  // Constructor for initialising all "combination trees"
  sharingtree(sharingtree<V> &S, sharingtree<V> &T)
      : maxLayers{std::max(S.maxLayers, T.maxLayers)},
        maxNodes{S.maxNodes + T.maxNodes}, maxSons{S.maxSons + T.maxSons} {
    layerBuffer = new st_layer[maxLayers];
    nodeBuffer = new st_node[maxNodes];
    sonBuffer = new st_son[maxSons];

    // Create root with arbitrary value, will not be read
    createNode(0, false, true);
    addFirstLayer();
  }

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
