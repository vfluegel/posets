#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include <posets/concepts.hh>
#include <posets/utils/sharingforest.hh>

namespace posets::downsets {

  template <Vector V>
  class simple_sharingtree_backed {
    private:
      size_t dim;
      size_t root {};
      std::shared_ptr<utils::sharingforest<V>> forest;
      std::vector<V> vector_set;
      static std::map<size_t, std::weak_ptr<utils::sharingforest<V>>> forest_map;

      void init_forest (size_t dimkey) {
        auto res = sharingtree_backed::forest_map.find (dimkey);
        if (res != sharingtree_backed::forest_map.end ()) {
          // two cases now, either the pointer is live and we take a lock on it
          // or it's dead, and we need to update it
          const std::shared_ptr<utils::sharingforest<V>> live = res->second.lock ();
          if (live) {
            this->forest = live;
          }
          else {
            this->forest = std::make_shared<utils::sharingforest<V>> (dimkey);
            res->second = this->forest;
          }
        }
        else {
          this->forest = std::make_shared<utils::sharingforest<V>> (dimkey);
          sharingtree_backed::forest_map.emplace (dimkey, this->forest);
        }
      }

      // Code borrowed from kdtree_backed, same idea as there to avoid
      // duplicates, and then to keep the antichain of max elements only
      void reset_tree (std::vector<V>&& elements) noexcept {
        std::vector<V*> pelements;
        pelements.reserve (elements.size ());
        for (auto& e : elements)
          pelements.push_back (&e);

        std::sort (pelements.begin (), pelements.end (),
                   // A strict total order.
                   [] (const V* v1, const V* v2) {
                     // A quite costly thing to do.
                     for (size_t i = 0; i < v1->size (); ++i) {
                       if ((*v1)[i] > (*v2)[i])
                         return false;
                       if ((*v1)[i] < (*v2)[i])
                         return true;
                     }
                     // Equal MUST return false.
                     return false;
                   });

        // then remove duplicates
        size_t dups_pos = pelements.size ();
        for (size_t i = pelements.size () - 1; i > 0; --i)
          if (*pelements[i] == *pelements[i - 1])
            std::swap (pelements[i], pelements[--dups_pos]);
        if (dups_pos != pelements.size ())
          pelements.erase (pelements.begin () + dups_pos, pelements.end ());

        // now, we can make a tree out of the set to eliminate dominated
        // elements
        std::vector<V> set_elements;
        set_elements.reserve (pelements.size ());
        for (auto& e : pelements)
          set_elements.push_back (e->copy ());

        size_t temp_tree = this->forest->add_vectors (std::move (set_elements), false);

        std::vector<V> antichain;
        antichain.reserve (pelements.size ());
        this->vector_set.clear ();
        this->vector_set.reserve (pelements.size ());
        for (auto& e : pelements) {
          if (not this->forest->covers_vector (temp_tree, *e, true)) {
            this->vector_set.push_back (e->copy ());
            antichain.push_back (std::move (*e));
          }
        }

        this->root = this->forest->add_vectors (std::move (antichain), false);
      }

    public:
      using value_type = V;

      sharingtree_backed () = delete;
      sharingtree_backed (const sharingtree_backed&) = delete;
      sharingtree_backed (sharingtree_backed&&) = default;
      sharingtree_backed& operator= (const sharingtree_backed&) = delete;
      sharingtree_backed& operator= (sharingtree_backed&&) = default;

      sharingtree_backed (std::vector<V>&& elements) {
        init_forest (elements.begin ()->size ());
        reset_tree (std::move (elements));
      }

      sharingtree_backed (V&& v) {
        init_forest (v.size ());
        this->root = this->forest->add_vectors (std::array<V, 1> {std::move (v)}, false);
        this->vector_set = this->forest->get_all (this->root);
      }

      [[nodiscard]] auto size () const { return this->vector_set.size (); }
      auto begin () { return this->vector_set.begin (); }
      [[nodiscard]] auto begin () const { return this->vector_set.begin (); }
      auto end () { return this->vector_set.end (); }
      [[nodiscard]] auto end () const { return this->vector_set.end (); }
      [[nodiscard]] auto& get_backing_vector () { return vector_set; }
      [[nodiscard]] const auto& get_backing_vector () const { return vector_set; }

      [[nodiscard]] bool contains (const V& v) const {
        return this->forest->covers_vector (this->root, v);
      }

      // Union in place
      void union_with (sharingtree_backed&& other) {
        const size_t op1 = this->root;
        const size_t op2 = other.root;
        const size_t new_root = this->forest->st_union (op1, op2);
        this->root = new_root;
        this->vector_set = this->forest->get_all (this->root);
        this->trim_by_dom ();

        assert (other.size () > 0);
        std::vector<V*> undomd;
        undomd.reserve (this->size () + other.size ());
        // for all elements in this tree, if they are not strictly
        // dominated by the other tree, we keep them
        for (auto& e : this->vector_set)
          if (not other.forest->covers_vector (other.root, e, true))
            undomd.push_back (&e);

        // for all elements in the other tree, if they are not dominated
        // (not necessarily strict) by this tree, we keep them
        for (auto& e : other.vector_set)
          if (not this->forest->covers_vector (this->root, e, false))
            undomd.push_back (&e);
        
        // ready to rebuild the tree now
        assert (not undomd.empty ());

        // We can move the referenced elements
        std::vector<V> result;
        result.reserve (undomd.size ());
        this->vector_set.clear ();
        this->vector_set.reserve (undomd.size ());
        for (const auto& r : undomd) {
          this->vector_set.push_back (r->copy ());
          result.push_back (std::move (*r));
        }

        this->root = this->forest->add_vectors (std::move (result), false);
      }

      // Intersection in place
      void intersect_with (const sharingtree_backed& other) {
        std::vector<V> intersection;
        bool smaller_set = false;

        for (const auto& x : this->vector_set) {
          assert (x.size () > 0);

          // If x is part of the set of all meets, then x will dominate the
          // whole list! So we use this to short-circuit the computation: we
          // first check whether x will be there (which happens only if it is
          // itself dominated by some element in other)
          const bool dominated = other.contains (x);
          if (dominated)
            intersection.push_back (x.copy ());
          else
            for (auto& y : other)
              intersection.push_back (x.meet (y));

          // If x wasn't in the set of meets, dominated is false and
          // the set of minima is different than what is in this->tree
          smaller_set or_eq not dominated;
        }

        // We can skip working all if we already have the antichain
        // of minimal elements
        if (not smaller_set)
          return;

        // Worst-case scenario: we do need to work
        this->reset_tree (std::move (intersection));
      }

      template <typename F>
      auto apply (const F& lambda) const {
        std::vector<V> ss;
        ss.reserve (this->vector_set.size ());

        for (const auto& v : this->vector_set)
          ss.push_back (lambda (v));

        return sharingtree_backed (std::move (ss));
      }
  };

  template <Vector V>
  std::map<size_t, std::weak_ptr<utils::sharingforest<V>>> sharingtree_backed<V>::forest_map;

  template <Vector V>
  inline std::ostream& operator<< (std::ostream& os, const sharingtree_backed<V>& f) {
    for (auto&& el : f.vector_set)
      os << el << std::endl;
    return os;
  }
}
