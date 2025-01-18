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
  class sharingtree_backed {
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

      void trim_by_dom () {
#ifdef SHARINGTREE_DOMTRIM
        std::vector<V*> refs;
        std::vector<V> undomd;
        refs.reserve (this->vector_set.size ());
        undomd.reserve (this->vector_set.size ());
        for (auto& v : this->vector_set)
          if (not this->forest->covers_vector (this->root, v, true))
            refs.push_back (&v);
        for (const auto& r : refs)
          undomd.push_back (std::move (*r));

# ifndef NDEBUG
        std::cout << "Undominated " << undomd.size () << "/" << this->vector_set.size () << '\n';
# endif

        this->root = this->forest->add_vectors (std::move (undomd));
        this->vector_set = this->forest->get_all (this->root);
#endif
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
        this->root = this->forest->add_vectors (std::move (elements));
        this->vector_set = this->forest->get_all (this->root);
        this->trim_by_dom ();
      }

      sharingtree_backed (V&& v) {
        init_forest (v.size ());
        this->root = this->forest->add_vectors (std::array<V, 1> {std::move (v)});
        this->vector_set = this->forest->get_all (this->root);
        this->trim_by_dom ();
      }

      [[nodiscard]] auto size () const { return this->vector_set.size (); }
      auto begin () { return this->vector_set.begin (); }
      [[nodiscard]] auto begin () const { return this->vector_set.begin (); }
      auto end () { return this->vector_set.end (); }
      [[nodiscard]] auto end () const { return this->vector_set.end (); }
      [[nodiscard]] auto& get_backing_vector () { return vector_set; }
      [[nodiscard]] const auto& get_backing_vector () const { return vector_set; }

      [[nodiscard]] bool contains (const V& v) const {
#ifdef SHARINGTREE_SIMCONTAINS
        size_t sappling = this->forest->add_vectors (std::array<V, 1> {v.copy ()});
        this->forest->print_children (sappling, 0);
        this->forest->print_children (this->root, 0);
        return this->forest->check_simulation (this->root, sappling);
#else
        return this->forest->covers_vector (this->root, v);
#endif
      }

      // Union in place
      void union_with (sharingtree_backed&& other) {
        const size_t op1 = this->root;
        const size_t op2 = other.root;
        const size_t new_root = this->forest->st_union (op1, op2);
        this->root = new_root;
        this->vector_set = this->forest->get_all (this->root);
        this->trim_by_dom ();
      }

      // Intersection in place
      void intersect_with (const sharingtree_backed& other) {
#ifdef SHARINGTREE_GRAPHINTER
        // Worst-case scenario: we do need to work
        this->root = this->forest->st_intersect (this->root, other.root);
        this->vector_set = this->forest->get_all (this->root);
#else
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
        this->root = this->forest->add_vectors (std::move (intersection));
        this->vector_set = this->forest->get_all (this->root);
#endif
        this->trim_by_dom ();
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
