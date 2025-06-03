#pragma once

#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <ranges>
#include <stack>
#include <tuple>
#include <vector>

#include <posets/concepts.hh>

namespace posets::utils {

// We will be (micro)managing (C-style) dynamic memory for the set of nodes we
// keep in each layer. The initial number of nodes allowed is determined by
// the number below (multiplied by the number of layers = the dimension + 1).
// When the buffer is full, we double its capacity and copy everything to the
// new block of reserved memory.
#ifndef SHARINGFOREST_INIT_LAYER_SIZE
# define SHARINGFOREST_INIT_LAYER_SIZE 100UL
#endif
#ifndef SHARINGFOREST_INIT_MAX_CHILDREN
# define SHARINGFOREST_INIT_MAX_CHILDREN 10UL
#endif

  // Forward definition for the operator<<
  template <Vector>
  class radix_sharingforest;

  template <Vector V>
  std::ostream& operator<< (std::ostream& os, const utils::radix_sharingforest<V>& f);

  template <Vector V>
  class radix_sharingforest {
    private:
      template <Vector V2>
      friend std::ostream& operator<< (std::ostream& os, const utils::radix_sharingforest<V2>& f);

      size_t dim;

      struct st_node {
          size_t numchild;
          std::vector<typename V::value_type> label;
          size_t cbuffer_offset;
      };

      // NOTE: We will be keeping a unique table of st_nodes using an unordered
      // map, for this we need a hash function and an equivalence operation.
      // Since nodes keep offsets instead of pointers, we need both the hash and
      // equivalence operations to be able to access the unique table. Hence, we
      // keep a reference to the table in them.
      class st_hash {
          radix_sharingforest* f;

        public:
          st_hash (radix_sharingforest* that) : f {that} {}

          size_t operator() (const st_node& k) const {
            size_t res = std::hash<typename V::value_type> () (k.label[0]);
            for (size_t v = 1; v < k.label.size(); v++) {
              res ^= std::hash<typename V::value_type> () (k.label[v]) << (v + 1);
            }
            size_t* children = f->child_buffer + k.cbuffer_offset;
            for (size_t i = 0; i < k.numchild; i++)
              res ^= std::hash<size_t> () (children[i]) << (i + 1);
            return res;
          }
      };

      class st_equal {
          radix_sharingforest* f;

        public:
          st_equal (radix_sharingforest* that) : f {that} {}

          bool operator() (const st_node& lhs, const st_node& rhs) const {
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

      std::vector<std::vector<st_node>> layers;
      size_t* child_buffer;
      size_t cbuffer_size;
      size_t cbuffer_nxt;

      // This is a cache/hash-map to get in-layer node identifiers from their
      // signature
      std::vector<std::unordered_map<st_node, size_t, st_hash, st_equal>> inverse;
      // Second cache/hash-map to check for a pair of node(-identifiers) in a layer
      // whether there is a simulation relation in the left-to-right direction
      std::vector<std::unordered_map<std::pair<size_t, size_t>, bool,
                                     boost::hash<std::pair<size_t, size_t>>>>
          simulating;
      // One more for union
      std::vector<std::unordered_map<std::pair<size_t, size_t>, size_t,
                                     boost::hash<std::pair<size_t, size_t>>>>
          cached_union;

      void init (size_t dim) {
        this->dim = dim;
        layers.resize (dim + 1);

        for (size_t i = 0; i < dim + 1; i++) {
          inverse.emplace_back (SHARINGFOREST_INIT_LAYER_SIZE, st_hash (this), st_equal (this));
          simulating.emplace_back ();
          cached_union.emplace_back ();
        }

        cbuffer_size = SHARINGFOREST_INIT_LAYER_SIZE * SHARINGFOREST_INIT_MAX_CHILDREN;
        child_buffer = new size_t[cbuffer_size];
        cbuffer_nxt = 0;
      }

      std::optional<size_t> has_son (st_node& node, size_t child_layer, int val) {
        if (node.numchild == 0)
          return std::nullopt;

        int left = 0;
        int right = node.numchild - 1;
        size_t* children = child_buffer + node.cbuffer_offset;

        while (left <= right) {
          int mid = left + (right - left) / 2;
          assert (mid < static_cast<int>(node.numchild));
          typename V::value_type mid_val = layers[child_layer][children[mid]].label[0];

          if (mid_val == val)
            return mid;
          if (mid_val > val)
            left = mid + 1;
          else
            right = mid - 1;
        }

        return std::nullopt;
      }

      bool simulates (size_t n1idx, size_t n2idx, size_t layidx) {
        auto node_pair = std::make_pair (n1idx, n2idx);
        auto cached = simulating[layidx].find (node_pair);
        if (cached != simulating[layidx].end ())
          return cached->second;

        st_node n1 = layers[layidx][n1idx];
        st_node n2 = layers[layidx][n2idx];
        // If the node is in the last layer, we just check the labels
        if (layidx == this->dim)
          return n1.label[0] >= n2.label[0];
        size_t* n1_children = child_buffer + n1.cbuffer_offset;
        size_t* n2_children = child_buffer + n2.cbuffer_offset;

        auto largest_child_n1 = n1.label.size () > 1 ? n1.label[1] : layers[layidx + 1][n1_children[0]].label[0];
        auto smallest_child_n2 = n2.label.size () > 1 ? n2.label[1] : layers[layidx + 1][n2_children[n2.numchild - 1]].label[0];

        if (n1.label[0] < n2.label[0] or largest_child_n1 < smallest_child_n2) {
          // If the label of n1 is too small or its largest child is already smaller than n2's
          // smallest, we already know it can't simulate
          return false;
        }

        // Stack contains node and child ID of S, node and child ID of T, layer
        std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;
        current_stack.emplace (n1idx, 0, n2idx, 0, layidx);

        while (not current_stack.empty ()) {
          auto [n1idx, c1, n2idx, c2, layidx] = current_stack.top ();
#ifndef NDEBUG
          std::cout << "Simulation check: " << "Node " << layidx << "." << n1idx << " (c=" << c1
                    << ") vs. " << layidx << "." << n2idx << " (c=" << c2 << ")\n";
#endif
          current_stack.pop ();
          n1 = layers[layidx][n1idx];
          n2 = layers[layidx][n2idx];

          if (n1.label.size() > 1) {
            std::vector<size_t> candidates {n2idx};
            size_t current_layidx = layidx;
            bool all_simulated = true;
            // Compare every branch in n2 -> every branch needs to be smaller 
            for (auto d = n1.label.begin (); d != n1.label.end () - 1; ++d) {
              std::vector<size_t> new_candidates;
              for (size_t otheridx : candidates) {
                st_node other_node = layers[current_layidx][otheridx];
                if (other_node.label.size () > 1) {
                  bool simulates = std::equal (d, n1.label.end (), other_node.label.begin (), other_node.label.end (),
                                              [](auto val1, auto val2) { return val1 >= val2; });
                  if (not simulates) {
                    all_simulated = false;
                    break;
                  }
                }
                else {
                  size_t* other_children = child_buffer + other_node.cbuffer_offset;
                  // all children need to be smaller, so we only have to consider nodes where the greatest child is simulated
                  if (other_node.label[0] <= *d and layers[current_layidx + 1][other_children[0]].label[0] <= *(d + 1)) {
                    for (size_t c = 0; c < other_node.numchild; c++) {
                      new_candidates.push_back (other_children[c]);
                    }
                  }
                  else {
                    all_simulated = false;
                    break;
                  }
                }
              }
              if (not all_simulated) break;
              current_layidx++;
              candidates = new_candidates;
            }
            simulating[layidx][std::make_pair (n1idx, n2idx)] = all_simulated;
            if (not current_stack.empty ()) {
              auto [m1idx, d1, m2idx, d2, ell] = current_stack.top ();
              current_stack.pop ();
              if (all_simulated)
                // The node simulates the other - this branch is done
                current_stack.emplace (m1idx, 0, m2idx, d2 + 1, ell);
              else
                // We have to check another child
                current_stack.emplace (m1idx, d1 + 1, m2idx, d2, ell);
            }
          }
          else if (n2.label.size () > 1) {
            // search through the subtree n1 to see if one path simulates it
            // Stack with node ID, child ID and layer offset
            std::stack<std::tuple<size_t, size_t, size_t>> df_stack;
            assert(c1 == 0); // We shouldn't have checked children already
            df_stack.emplace (n1idx, 0, 0);
            bool found_simulating = false;
            while (not df_stack.empty ()) {
              auto [other_id, child, layer_offset] = df_stack.top ();
              df_stack.pop ();
              st_node other_node = layers[layidx + layer_offset][other_id];
              if (other_node.label.size () > 1) {
                bool simulates = std::equal (n2.label.begin () + layer_offset, n2.label.end (), other_node.label.begin (), other_node.label.end (),
                                              [](auto val2, auto val1) { return val1 >= val2; });
                if (simulates) {
                  found_simulating = true;
                  break;
                }
                else if (not df_stack.empty ()) {
                  // No hit yet... Investigate the next child
                  auto [prev_id, child, prev_layer] = df_stack.top ();
                  df_stack.pop ();
                  df_stack.emplace (prev_id, child + 1, prev_layer);
                }
              }
              else if (layidx + layer_offset == this->dim) {
                // We have reached the end - the node is simulated
                found_simulating = true;
                break;
              }
              else if (child == other_node.numchild and not df_stack.empty ()) {
                // All children were checked but no luck - move on to next subtree
                auto [prev_id, child, prev_layer] = df_stack.top ();
                df_stack.pop ();
                df_stack.emplace (prev_id, child + 1, prev_layer);
              }
              else if (child < other_node.numchild) {
                size_t* other_children = child_buffer + other_node.cbuffer_offset;
                if (layers[layidx + layer_offset + 1][other_children[child]].label[0] >= n2.label[layer_offset + 1]) {
                  // This is a candidate - we have to go down this subtree
                  df_stack.emplace (other_id, child, layer_offset);
                  df_stack.emplace (other_children[child], 0, layer_offset + 1);
                }
                else {
                  // We move on to the next
                  df_stack.emplace (other_id, child + 1, layer_offset);
                }

              }
            }
            simulating[layidx][std::make_pair (n1idx, n2idx)] = found_simulating;
            if (not current_stack.empty ()) {
              auto [m1idx, d1, m2idx, d2, ell] = current_stack.top ();
              current_stack.pop ();
              if (found_simulating)
                // The node simulates the other - this branch is done
                current_stack.emplace (m1idx, 0, m2idx, d2 + 1, ell);
              else
                // We have to check another child
                current_stack.emplace (m1idx, d1 + 1, m2idx, d2, ell);
            }
          }
          // This is our base case, no more things to check on the n2 side
          // We now claim success for n2 being simulated by n1 and update the top
          // of the stack so that when we go back up in the tree we have one less
          // branch to check on the n2 side
          else if ((c2 == n2.numchild and n2.label.size() == 1) or layidx == this->dim) {
            assert (c2 == n2.numchild);
            simulating[layidx][std::make_pair (n1idx, n2idx)] = true;
            if (not current_stack.empty ()) {
              auto [m1idx, d1, m2idx, d2, ell] = current_stack.top ();
              current_stack.pop ();
              current_stack.emplace (m1idx, 0, m2idx, d2 + 1, ell);
            }
          }
          // Another base case: we've iterated through all the children on the
          // n1 side and failed to find a simulating one
          else if (c1 == n1.numchild and n1.label.size() == 1) {
#ifndef NDEBUG
            std::cout << "Another base case (layidx=" << layidx << ")\n";
#endif
            simulating[layidx][std::make_pair (n1idx, n2idx)] = false;
            if (not current_stack.empty ()) {
              auto [m1idx, d1, m2idx, d2, ell] = current_stack.top ();
              current_stack.pop ();
              current_stack.emplace (m1idx, d1 + 1, m2idx, d2, ell);
            }
          }
          // The last case is that we have two valid children indices,
          // then we have to go deeper in the product tree (lest the cache saves
          // us, or we find that the subtree will just certainly not satisfy
          // simulation)
          else {
            n1_children = child_buffer + n1.cbuffer_offset;
            n2_children = child_buffer + n2.cbuffer_offset;
            node_pair = std::make_pair (n1_children[c1], n2_children[c2]);
            cached = simulating[layidx + 1].find (node_pair);
            // Did we get lucky with the cache? then push back an updated node
            // with less obligations or keep searching on the n1 side
            if (cached != simulating[layidx + 1].end ()) {
#ifndef NDEBUG
              std::cout << "Got lucky with simulate cache!\n";
#endif
              if (cached->second)
                current_stack.emplace (n1idx, 0, n2idx, c2 + 1, layidx);
              else
                current_stack.emplace (n1idx, c1 + 1, n2idx, c2, layidx);
            }
            else {
              const st_node c1_node = layers[layidx + 1][n1_children[c1]];
              const st_node c2_node = layers[layidx + 1][n2_children[c2]];
              size_t* c1_children = child_buffer + c1_node.cbuffer_offset;
              size_t* c2_children = child_buffer + c2_node.cbuffer_offset;
              // In some cases, we can eliminate the subtree altogether
              bool children_disqualify = false;
              if (layidx + 1 < this->dim) {
                auto largest_child_c1 = c1_node.label.size() > 1 ? c1_node.label[1] : layers[layidx + 2][c1_children[0]].label[0];
                auto smallest_child_c2 = c2_node.label.size() > 1 ? c2_node.label[1] : layers[layidx + 2][c2_children[c2_node.numchild - 1]].label[0];
                children_disqualify = largest_child_c1 < smallest_child_c2;
              }
              
              if (c1_node.label[0] < c2_node.label[0] or children_disqualify) {
                current_stack.emplace (n1idx, c1 + 1, n2idx, c2, layidx);
                // Alright, we have to go deeper now; and push back what we just
                // popped
              }
              else {
                current_stack.emplace (n1idx, c1, n2idx, c2, layidx);
                current_stack.emplace (n1_children[c1], 0, n2_children[c2], 0, layidx + 1);
              }
            }
          }
        }
        node_pair = std::make_pair (n1idx, n2idx);
        cached = simulating[layidx].find (node_pair);
        assert (cached != simulating[layidx].end ());
        return cached->second;
      }

      void add_son (st_node& node, size_t son_layer, size_t son) {
        // Only add children if the label is a single item long
        assert (node.label.size() == 1);
        const int last = node.numchild - 1;
        const st_node& son_node = layers[son_layer][son];
        size_t* children = child_buffer + node.cbuffer_offset;

        // the new node is smaller than all existing ones
        assert (last == -1 or layers[son_layer][children[last]].label[0] > son_node.label[0]);

        // Insert the new value
        children[last + 1] = son;
        node.numchild++;
      }

      void add_son_unordered (st_node& node, size_t son_layer, size_t son) {
        assert (son_layer <= this->dim);  // Can't add son after the last layer
        // Find the insertion point using binary search
        int left = 0;
        int right = node.numchild - 1;
        const st_node& son_node = layers[son_layer][son];
        size_t* children = child_buffer + node.cbuffer_offset;
        while (left <= right) {
          const int mid = left + (right - left) / 2;
          assert (mid < static_cast<int> (node.numchild));
          assert (mid >= 0);
          assert (children[mid] < layers[son_layer].size ());
          const typename V::value_type mid_val = layers[son_layer][children[mid]].label[0];

          if (son_node.label[0] == mid_val) {
            const size_t new_son = node_union (son, children[mid], son_layer);
            children = child_buffer + node.cbuffer_offset;
            children[mid] = new_son;
            return;
          }
          if (mid_val < son_node.label[0])
            right = mid - 1;
          else
            left = mid + 1;
        }
        // Shift elements in the child buffer to make room for the new child
        for (int i = node.numchild; i > left; i--)
          children[i] = children[i - 1];

        // Insert the new value
        children[left] = son;
        node.numchild++;
      }

      void add_son_if_not_simulated (size_t node, size_t son_layer, st_node& father) {
        size_t* siblings = child_buffer + father.cbuffer_offset;
        for (size_t s = 0; s < father.numchild; s++)
          if (simulates (siblings[s], node, son_layer)) {
            return;
          }
            
        add_son (father, son_layer, node);
      }

      size_t add_node (st_node& node, size_t destination_layer) {
        auto existing_node = inverse[destination_layer].find (node);
        if (existing_node == inverse[destination_layer].end ()) {
          const size_t new_id = layers[destination_layer].size ();
          inverse[destination_layer][node] = new_id;
          layers[destination_layer].push_back (node);
          return new_id;
        }
        return existing_node->second;
      }

      size_t add_children (int num_child) {
        // Double the child buffer if it is full
        if (cbuffer_nxt + num_child >= cbuffer_size) {
          auto* new_buffer = new size_t[cbuffer_size * 2];
          std::memcpy (new_buffer, child_buffer, cbuffer_size * sizeof (size_t));
          cbuffer_size *= 2;
          delete[] child_buffer;
          child_buffer = new_buffer;
        }
        const size_t res = cbuffer_nxt;
        cbuffer_nxt += num_child;
        return res;
      }

      bool is_simulated (size_t nodeidx, st_node& father, size_t destination_layer) {
        size_t* siblings = child_buffer + father.cbuffer_offset;
        for (size_t s = 0; s < father.numchild; s++)
          if (simulates (siblings[s], nodeidx, destination_layer))
            return true;
        return false;
      }

      std::pair<size_t, bool> add_if_not_simulated (st_node& node, size_t destination_layer,
                                                    st_node& father) {
        const size_t res = add_node (node, destination_layer);
        return std::make_pair (res, is_simulated (res, father, destination_layer));
      }

      size_t node_union (size_t ns, size_t nt, size_t destination_layer) {
        // Stack contains node and child ID of S, node and child ID of T, layer
        std::stack<std::tuple<size_t, size_t, size_t, size_t, size_t>> current_stack;

        const st_node root_nod_e1 = layers[destination_layer][ns];
        const st_node root_nod_e2 = layers[destination_layer][nt];
        assert (root_nod_e1.numchild > 0 and root_nod_e2.numchild > 0);

        current_stack.emplace (ns, 0, nt, 0, destination_layer);

        while (not current_stack.empty ()) {
          auto [n_s, c_s, n_t, c_t, layer] = current_stack.top ();
          current_stack.pop ();
          assert (n_s < layers[layer].size ());
          const st_node node_s = layers[layer][n_s];
          assert (n_t < layers[layer].size ());
          const st_node node_t = layers[layer][n_t];

          // It's the first time we see this combination of nodes, so let's insert a draft
          // of it, still dirty, will be cleaned later
          if (c_s == 0 and c_t == 0) {
            assert (node_s.label[0] == node_t.label[0]);
            // Before creating a draft node and continuing, let's check the
            // cache
            auto cache_res = cached_union[layer].find (std::make_pair (n_s, n_t));
            if (layer > destination_layer and cache_res != cached_union[layer].end ()) {
#ifndef NDEBUG
              std::cout << "Avoided node in union = cache hit.\n";
#endif
              auto& father = layers[layer - 1].back ();
              if (not is_simulated (cache_res->second, father, layer))
                add_son (father, layer, cache_res->second);
              continue;
            }
            std::vector<typename V::value_type> new_label {node_s.label[0]};
            // Not found, so draft a node up
            if (layer < this->dim) {
              size_t n_s_children = node_s.label.size() > 1 ? 2 : node_s.numchild;
              size_t n_t_children = node_t.label.size() > 1 ? 2 : node_t.numchild;
              layers[layer].emplace_back (0, new_label, 
                                          add_children (n_s_children + n_t_children));
            }
            else {
              layers[layer].emplace_back (0, new_label);
            }
          }

          auto& under_construction = layers[layer].back ();
          if (node_s.label.size () > 1 or node_t.label.size () > 1) {
            assert(c_s == 0 and c_t == 0);
            auto& compressed_node_label = node_s.label.size () > 1 ? node_s.label : node_t.label;
            auto other_node_id = node_s.label.size () > 1 ? n_t : n_s;
            auto& other_node = layers[layer][other_node_id];
            
            std::vector<size_t> equal_path {other_node_id};
            std::optional<size_t> identical_son {};
            if (other_node.numchild > 0) {
              identical_son = has_son (other_node, layer + 1, compressed_node_label[1]);
            }
            // First go down the tree: As long as we find equal children, we add them to the path
            while (identical_son.has_value ()) {
              std::vector<typename V::value_type> identical_label {compressed_node_label[equal_path.size ()]};
              size_t prev_identical = identical_son.value ();
              auto identical_node = layers[layer + equal_path.size ()][prev_identical];
              layers[layer + equal_path.size ()].emplace_back (0, identical_label, add_children (identical_node.numchild + 1));
              
              equal_path.push_back (prev_identical);
              if (identical_node.numchild > 0 and layer + equal_path.size () < this->dim) {
                identical_son = has_son (identical_node, layer + equal_path.size (), compressed_node_label[equal_path.size ()]);
              }
              else 
                identical_son.reset ();
            }

            assert (layer + equal_path.size () <= this->dim);
            size_t ndix_to_union = equal_path.back ();
            st_node& node_to_union = layers[layer + equal_path.size () - 1][ndix_to_union];
            st_node union_son = layers[layer + equal_path.size () - 1].back ();

            if (layer + equal_path.size() == this->dim) {
              // The compressed node is completely equal to a path in the other node - we just use that
              auto& father = layers[layer - 1].back ();
              layers[layer].pop_back ();
              
              add_son (father, layer, other_node_id);
              cached_union[layer][std::make_pair (n_s, n_t)] = other_node_id;
              continue;
            }

            // Node still good here!!
            if (node_to_union.label.size () > 1) {
              // Special case: we have a node with a multi-element label (can only happen for the last!)
              bool compressed_simulates = std::equal (compressed_node_label.begin () + equal_path.size (), compressed_node_label.end (), node_to_union.label.begin () + 1, node_to_union.label.end (),
                                              [](auto compressed, auto other) { return compressed >= other; });
              if (compressed_simulates) {
                // We only need to keep the original vector
                union_son.label.insert (union_son.label.end (), compressed_node_label.begin () + equal_path.size (), compressed_node_label.end ());
              }
              else {
                // Also check the other direction
                bool other_simulates = std::equal (compressed_node_label.begin () + equal_path.size (), compressed_node_label.end (), node_to_union.label.begin () + 1, node_to_union.label.end (),
                                          [](auto compressed, auto other) { return compressed <= other; });
                if (other_simulates) {
                  union_son.label.insert (union_son.label.end (), node_to_union.label.begin () + 1, node_to_union.label.end ());
                }
                else {
                  // The two multi-element vectors might have more equivalent elements, we have to create those
                  std::vector<st_node> intermediaries {};
                  size_t i;
                  size_t length_offset = equal_path.size ();
                  assert (node_to_union.label.size () > 1);
                  assert (compressed_node_label.size () > length_offset);
                  for (i = 1; compressed_node_label[length_offset] == node_to_union.label[i]; i++) {
                    // The two vectors are equal - we need intermediary nodes
                    std::vector<typename V::value_type> intermediary_label {compressed_node_label[length_offset]};
                    intermediaries.emplace_back (0, intermediary_label, add_children (2));
                    length_offset++;
                    assert (i + 1 < node_to_union.label.size ());
                    assert (length_offset < compressed_node_label.size ());
                  }
                  assert (layer + length_offset <= this->dim);
                  // Create the two new sons
                  std::vector<typename V::value_type> new_compressed_label {compressed_node_label.begin () + length_offset, compressed_node_label.end()};
                  st_node new_compressed_son {0, new_compressed_label};
                  size_t compressed_son_id = add_node (new_compressed_son, layer + length_offset);
                  std::vector<typename V::value_type> new_other_label {node_to_union.label.begin () + i, node_to_union.label.end ()};
                  st_node new_other_son {0, new_other_label};
                  size_t other_son_id = add_node (new_other_son, layer + length_offset);
                  // Add the nodes as son (in order)
                  if (intermediaries.size () > 0) {
                    st_node intermediary_father = intermediaries.back ();
                    intermediaries.pop_back ();
                    if (new_compressed_label[0] > new_other_label[0]) {
                      add_son (intermediary_father, layer + length_offset, compressed_son_id);
                      add_son (intermediary_father, layer + length_offset, other_son_id);
                    }
                    else {
                      add_son (intermediary_father, layer + length_offset, other_son_id);
                      add_son (intermediary_father, layer + length_offset, compressed_son_id);
                    }
                    length_offset--;
                    size_t last_id = add_node (intermediary_father, layer + length_offset);
                    while (intermediaries.size () > 0) {
                      add_son (intermediaries.back (), layer + length_offset, last_id);
                      length_offset--;
                      last_id = add_node (intermediaries.back (), layer + length_offset);
                      intermediaries.pop_back ();
                    }
                    add_son(union_son, layer + length_offset, last_id);
                  }
                  else {
                    if (new_compressed_label[0] > new_other_label[0]) {
                      add_son (union_son, layer + length_offset, compressed_son_id);
                      add_son (union_son, layer + length_offset, other_son_id);
                    }
                    else {
                      add_son (union_son, layer + length_offset, other_son_id);
                      add_son (union_son, layer + length_offset, compressed_son_id);
                    }
                  }
                }
              }
            }
            else {
              // Create the new node with the remaining label
              std::vector<typename V::value_type> remaining_compressed_label {compressed_node_label.begin () + equal_path.size (), 
                                                                            compressed_node_label.end ()};
              st_node new_compressed_node {0, remaining_compressed_label};
              size_t new_node_idx = add_node (new_compressed_node, layer + equal_path.size ());

              // Add the sons to the union_son node
              size_t child_layer = layer + equal_path.size ();
              size_t current_child = 0;
              size_t* original_children = child_buffer + node_to_union.cbuffer_offset;
              if (node_to_union.numchild > 0) {
                st_node& original_child = layers[child_layer][original_children[0]];
                // Add all the children with larger labels first
                while (current_child < node_to_union.numchild and original_child.label[0] > new_compressed_node.label[0]) {
                  add_son (union_son, child_layer, original_children[current_child]);
                  original_child = layers[child_layer][original_children[current_child]];
                  current_child++;
                }
              }
              add_son_if_not_simulated (new_node_idx, child_layer, union_son);
              for (; current_child < node_to_union.numchild; current_child++) {
                add_son_if_not_simulated (original_children[current_child], child_layer, union_son);
              }
            }

            auto& father_path = layers[layer - 1].back ();
            size_t* father_pathchildren = father_path.cbuffer_offset + child_buffer;
            for (size_t c = 0; c < father_path.numchild; c++) {
              print_children (father_pathchildren[c], layer);
            }

            layers[layer + equal_path.size () - 1].pop_back ();
            
            size_t union_son_idx = add_node (union_son, layer + equal_path.size () - 1);
            // Now its time to go back up 
            while (equal_path.size () > 1) {
              equal_path.pop_back ();
              size_t ndix_to_union = equal_path.back ();
              st_node& node_to_union = layers[layer + equal_path.size () - 1][ndix_to_union];
              
              st_node equal_father = layers[layer + equal_path.size () - 1].back ();
              layers[layer + equal_path.size () - 1].pop_back ();
              size_t child_layer = layer + equal_path.size ();
              size_t current_child = 0;
              size_t* original_children = child_buffer + node_to_union.cbuffer_offset;
              if (node_to_union.numchild > 0) {
                st_node& original_child = layers[child_layer][original_children[0]];
                // Add all the children with larger labels first
                while (current_child < node_to_union.numchild and original_child.label[0] > union_son.label[0]) {
                  add_son (equal_father, child_layer, original_children[current_child]);
                  original_child = layers[child_layer][original_children[current_child]];
                  current_child++;
                }
              }
              // Skip the equal child - we replace it with our new one
              current_child++;
              // Then add the new son
              add_son_if_not_simulated (union_son_idx, child_layer, equal_father);
              // And also add the smaller children
              for (; current_child < node_to_union.numchild; current_child++) {
                add_son_if_not_simulated (original_children[current_child], child_layer, equal_father);
              }
              union_son_idx = add_node (equal_father, child_layer - 1);
              union_son = layers[child_layer - 1][union_son_idx];
            }

            assert(layer + equal_path.size () - 1 == layer);
            st_node created = layers[layer][union_son_idx];
            // Lastly, we commit the node
            auto& father = layers[layer - 1].back ();
            cached_union[layer][std::make_pair (n_s, n_t)] = union_son_idx;
            add_son (father, layer, union_son_idx);

          }
          // Base case: We are done with this combo and can clean up
          else if (c_s == node_s.numchild and c_t == node_t.numchild) {
            if (layer == destination_layer) {
              // We skip this step, we clean the root after the loop, but we are done
              break;
            }
            // The very first thing to do is to check whether our draft of node is
            // not a repetition of something in the table already!
            auto& father = layers[layer - 1].back ();
            layers[layer].pop_back ();
            auto [union_res, domd] = add_if_not_simulated (under_construction, layer, father);
            cached_union[layer][std::make_pair (n_s, n_t)] = union_res;
            if (not domd)
              add_son (father, layer, union_res);
            // Recursive step: Either just add the son to the draft node and
            // continue with the next node, or continue down the tree if
            // necessary
          }
          else if (layer < this->dim) {
            size_t* node_s_children = child_buffer + node_s.cbuffer_offset;
            size_t* node_t_children = child_buffer + node_t.cbuffer_offset;

            // Case 1: One node is "done"
            // We add the existing node (including all its sons!) to the node
            // currently being constructed
            if (node_s.label.size () == 1 and c_s == node_s.numchild) {
              add_son_if_not_simulated (node_t_children[c_t], layer + 1, under_construction);
              if (c_t < node_t.numchild)
                current_stack.emplace (n_s, c_s, n_t, c_t + 1, layer);
            }
            else if (node_t.label.size () == 1 and c_t == node_t.numchild) {
              add_son_if_not_simulated (node_s_children[c_s], layer + 1, under_construction);
              if (c_s < node_s.numchild)
                current_stack.emplace (n_s, c_s + 1, n_t, c_t, layer);
            }
            else {
              const st_node& son_s = layers[layer + 1][node_s_children[c_s]];
              const st_node& son_t = layers[layer + 1][node_t_children[c_t]];
              // Case 2: One child is "ahead", We add just as when one list is done
              if (son_t.label[0] > son_s.label[0]) {
                add_son_if_not_simulated (node_t_children[c_t], layer + 1, under_construction);
                if (c_t < node_t.numchild)
                  current_stack.emplace (n_s, c_s, n_t, c_t + 1, layer);
              }
              else if (son_s.label[0] > son_t.label[0]) {
                add_son_if_not_simulated (node_s_children[c_s], layer + 1, under_construction);
                if (c_s < node_s.numchild)
                  current_stack.emplace (n_s, c_s + 1, n_t, c_t, layer);
              }
              // Case 3: The values of the sons match, so we need to continue the recursion
              else {
                assert (son_s.label[0] == son_t.label[0]);
                if (c_s < node_s.numchild and c_t < node_t.numchild)
                  current_stack.emplace (n_s, c_s + 1, n_t, c_t + 1, layer);
                current_stack.emplace (node_s_children[c_s], 0, node_t_children[c_t], 0,
                                       layer + 1);
              }
            }
          }
        }

        // Clean up the root, we remove it and re-add it, so it is checked whether an identical
        // node exists
        auto constructed_node = layers[destination_layer].back ();
        layers[destination_layer].pop_back ();
        return add_node (constructed_node, destination_layer);
      }

      // NOLINTBEGIN(misc-no-recursion)
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
      size_t build_node (std::vector<size_t>& vecs, size_t current_layer,
                         const auto& element_vec) {
        assert (not vecs.empty ());
        st_node new_node {0};
        if (current_layer > 0 and vecs.size () == 1) {
          auto& vec = element_vec[vecs[0]];
          std::vector<typename V::value_type> label_vec {vec.begin () + (current_layer - 1), vec.end ()};
          new_node.label = label_vec;
        }
        else {
          // If currentLayer is 0, we set the label to the dummy value -1 for the root
          // Else all nodes should have the same value at index currentLayer - 1, so
          // we just use the first
          const typename V::value_type label {current_layer == 0
                                                  ? static_cast<typename V::value_type> (-1)
                                                  : element_vec[vecs[0]][current_layer - 1]};
          
          std::vector<typename V::value_type> label_vec {label};
          new_node.label = label_vec;
          
          // We have not reached the last layer - so add children
          if (current_layer < this->dim) {
            // Partition and order the future children
            // TODO: Try to do constant-access bucketing based on the value of k.
            // Probably won't pay off unless the set of vectors we are adding is
            // dense in most components.
            std::map<typename V::value_type, std::vector<size_t>, std::greater<>> new_partition {};
            for (const auto& vec : vecs)
              new_partition[element_vec[vec][current_layer]].push_back (vec);
            new_node.cbuffer_offset = add_children (new_partition.size ());

            for (auto& [n, children] : new_partition) {
              // Build a new son for each individual value at currentLayer + 1
              const size_t new_son = build_node (children, current_layer + 1, element_vec);
              bool found = false;
              size_t* current_children = child_buffer + new_node.cbuffer_offset;
              for (size_t s = 0; s < new_node.numchild; s++) {
                if (simulates (current_children[s], new_son, current_layer + 1)) {
                  found = true;
                  break;
                }
              }
              if (not found)
                add_son (new_node, current_layer + 1, new_son);
            }
          }
        }
        return add_node (new_node, current_layer);
      }
      // NOLINTEND(misc-no-recursion)

    public:
      radix_sharingforest () = default;

      radix_sharingforest (size_t dim) { this->init (dim); }

      ~radix_sharingforest () {
        if (layers.empty ())
          return;

        delete[] child_buffer;
      }

      [[nodiscard]] std::vector<V> get_all (std::optional<size_t> root = {}) const {
        // Stack with tuples (layer, node id, child id)
        std::stack<std::tuple<size_t, size_t, size_t>> to_visit;

        if (root) {
          const st_node& root_node = layers[0][root.value ()];
          size_t* first_children = child_buffer + root_node.cbuffer_offset;
          for (size_t c = 0; c < root_node.numchild; c++)
            to_visit.emplace (1, first_children[c], 0);
        }
        else {
          // Add all roots at dimension 0
          for (size_t i = 0; i < layers[1].size (); i++)
            to_visit.emplace (1, i, 0);
        }
        
        std::vector<V> res;
        std::vector<typename V::value_type> temp;
        while (not to_visit.empty ()) {
          const auto [lay, node, child] = to_visit.top ();
          to_visit.pop ();
          const auto parent = layers[lay][node];
          // first time we see parent? get its label
          if (child == 0)
            temp.insert(temp.end (), parent.label.begin (), parent.label.end ());

          // base case: reached the bottom layer
          if (lay == this->dim or parent.label.size () > 1) {
            assert (child == 0);
            std::vector<typename V::value_type> cpy {temp};
            res.push_back (V (std::move (cpy)));
            temp.resize (temp.size () - parent.label.size ());  // done with this node
          }
          else {  // recursive case
            // Either we're done with this node and we just mark it as visited or
            // we need to keep it and we add it's next son
            size_t* children = child_buffer + parent.cbuffer_offset;
            if (child < parent.numchild) {
              to_visit.emplace (lay, node, child + 1);
              to_visit.emplace (lay + 1, children[child], 0);
            }
            else {
              temp.pop_back ();  // done with this parent
            }
          }
        }
        return res;
      }

      void print_children (size_t n, size_t layer) {
#ifndef NDEBUG
        assert (layer <= this->dim);
        st_node& node = layers[layer][n];
        size_t* children = child_buffer + node.cbuffer_offset;
        std::cout << std::string (layer, '\t') << layer << "." << n << " [";
        for (auto val : node.label)
          std::cout << static_cast<int> (val) << " ";
        std::cout << "] -> (" << (layer == this->dim ? "" : "\n");
        for (size_t i = 0; i < node.numchild; i++)
          print_children (children[i], layer + 1);
        std::cout << (layer == this->dim ? " " : std::string (layer, '\t')) << " )\n";
#endif
      }

      size_t st_union (size_t root1, size_t root2) { return node_union (root1, root2, 0); }

      /* Recursive domination check of given vector by vectors in the language of
       * the given tree root.
       *
       * Complexity: The implementation below has complexity O(n.d) where n is
       * the number of vectors, d is the number of dimensions. This is the same as
       * list-based data structures for vectors. However, since it is a DFS of the
       * DFA representation of (bisimulation non-dominated) vectors, it could be
       * exponentially faster than this.
       */
      bool covers_vector (size_t root, const V& covered, bool strict = false) {
        // Stack with tuples (layer, node id, strictness, child id)
        std::stack<std::tuple<size_t, size_t, bool, size_t>> to_visit;
        // Visited cache
        std::vector<std::unordered_map<size_t, bool>> visited;
        for (size_t i = 0; i < this->dim + 1; i++)
          visited.emplace_back ();

        // Add all roots at dimension 0 such that their labels cover the first
        // component of the given vector
        const st_node& root_node = layers[0][root];
        size_t* root_children = child_buffer + root_node.cbuffer_offset;
        for (size_t i = 0; i < root_node.numchild; i++) {
          assert (root_children[i] < layers[1].size ());
          if (layers[1][root_children[i]].label.size() > 1) {
            auto child_label = layers[1][root_children[i]].label;
            // The child already is condensed, we can compare the entire vector
            bool covers = std::equal(child_label.begin(), child_label.end(), covered.begin(), covered.end(),
                                              [strict](auto child_val, auto covered_val) { return strict ? (child_val > covered_val) : (child_val >= covered_val); });
            if (covers)
              return true;
          }
          else if (covered[0] <= layers[1][root_children[i]].label[0]) {
            const bool owe_strict = strict and covered[0] == layers[1][root_children[i]].label[0];
            to_visit.emplace (1, root_children[i], owe_strict, 0);
          }
        }

        while (not to_visit.empty ()) {
          const auto [lay, node, owe_strict, child] = to_visit.top ();
          assert (lay <= this->dim);
          assert (node < layers[lay].size ());
          to_visit.pop ();
          // if this node has been visited, it means this node is useless
          // NOTE: this works only because we are doing a DFS and not a BFS
          if (child == 0) {
            auto res = visited[lay].find (node);
            // if we visited with strictness then we can safely exit only if
            // we come back with strictness prev_strict -> owe_strict, i.e.
            // not prev_strict or owe_strict
            if (res != visited[lay].end ()) {
              if ((not res->second) or owe_strict) {
#ifndef NDEBUG
                std::cout << "Avoided node in covers check, DFS cache helps.\n";
#endif
                continue;
              }
              // we conjoin with previous result if any to make sure we have
              // more early exits based on implication condition above
              visited[lay][node] = res->second and owe_strict;
            }
            else {
              // no early exit? then mark the node as visited and keep going
              visited[lay][node] = owe_strict;
            }
          }
          const auto parent = layers[lay][node];

          // base case: reached the bottom layer
          if (lay == this->dim) {
            assert (child == 0);
            if (covered[lay - 1] < parent.label[0] or
                ((not owe_strict) and covered[lay - 1] == parent.label[0]))
              return true;
          }
          else if (parent.label.size () > 1) {
            // We found a condensed node - we can compare all the rest
            bool covers = std::equal (parent.label.begin (), parent.label.end (), covered.begin () + (lay - 1), covered.end (),
                                              [strict, owe_strict](auto node_val, auto covered_val) { return (strict && owe_strict) ? (node_val > covered_val) : (node_val >= covered_val); });
            if (covers)
              return true;
          }
          else {  // recursive case
            // Either we're done with this node and we just mark it as visited or
            // we need to keep it and we add it's next son
            assert (child < parent.numchild);
            size_t* children = child_buffer + parent.cbuffer_offset;
            const size_t c = child;
            auto child_node = layers[lay + 1][children[c]];
            // early exit if the largest child is smaller
            if (covered[lay] > child_node.label[0] or
                (owe_strict and covered[lay] >= child_node.label[0]))
              continue;

            const bool still_owe_strict = owe_strict and covered[lay] == child_node.label[0];
            assert (c < parent.numchild);
            if (c + 1 < parent.numchild)
              to_visit.emplace (lay, node, owe_strict, c + 1);
            to_visit.emplace (lay + 1, children[c], still_owe_strict, 0);
          }
        }
        return false;
      }

      template <std::ranges::input_range R>
      size_t add_vectors (R&& elements) {
        assert (not layers.empty ());

        auto element_vec = std::forward<R> (elements);
        // We start a Trie encoded as a map from prefixes to sets of indices of
        // the original vectors, we insert the root too: an empty prefix mapped to
        // the set of all indices
        std::vector<size_t> vector_ids (element_vec.size ());
        std::iota (vector_ids.begin (), vector_ids.end (), 0);

        const size_t root_id = build_node (vector_ids, 0, element_vec);
#ifndef NDEBUG
        size_t maxlayer = 0;
        size_t totlayer = 0;
        for (size_t i = 0; i < this->dim + 1; i++) {
          totlayer += this->layers[i].size ();
          if (maxlayer < this->layers[i].size ())
            maxlayer = this->layers[i].size ();
        }
        std::cout << "[" << this->dim << " Forest stats] Max layer size=" << maxlayer
                  << " total layers' size=" << totlayer << " (bytes per el=" << sizeof (st_node)
                  << ")" << '\n';
        std::cout << "[" << this->dim << " Forest stats] Child buff size=" << this->cbuffer_size
                  << " (bytes per el=" << sizeof (size_t) << ")\n";
#endif
        return root_id;
      }

      /*
        For testing: Check that all children are ordered descending
      */
      bool check_child_order () {
        size_t layer_num = 0;
        for (auto& l : layers) {
          for (st_node& n : l) {
            size_t* children = child_buffer + n.cbuffer_offset;
            for (size_t i = 1; i < n.numchild; i++) {
              // There is a child with a larger label than the previous child
              if (layers[layer_num + 1][children[i]].label[0] >
                  layers[layer_num + 1][children[i - 1]].label[0]) {
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
      bool check_simulation (size_t n1, size_t n2) { return simulates (n1, n2, 0); }
  };

  template <Vector V>
  inline std::ostream& operator<< (std::ostream& os, const radix_sharingforest<V>& f) {
    for (auto&& el : f.get_all ())
      os << el << '\n';

    os << "Layers:" << '\n';
    for (auto& l : f.layers) {
      for (auto& n : l) {
        os << "[ ";
        for (auto& v : n.label)
          os << static_cast<int> (v) << " ";
        os << "]";
      } 
      os << '\n';
    }
    return os;
  }

}  // namespace posets::utils
