#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <posets/concepts.hh>

namespace posets::downsets {
  // A forward definition to allow for friend status
  template <Vector V>
  class vector_or_kdtree_backed;

  template <Vector V>
  class vector_backed {
    public:
      using value_type = V;

      vector_backed (V&& v) { insert (std::move (v)); }

      vector_backed (std::vector<V>&& elements) noexcept {
        assert (not elements.empty ());
        for (auto&& e : elements)
          insert (std::move (e));
      }

    private:
      vector_backed () = default;
      std::vector<V> vector_set;

    public:
      vector_backed (const vector_backed&) = delete;
      vector_backed (vector_backed&&) = default;
      vector_backed& operator= (vector_backed&&) = default;
      vector_backed& operator= (const vector_backed&) = delete;

      bool operator== (const vector_backed& other) = delete;

      [[nodiscard]] bool contains (const V& v) const {
        return std::ranges::any_of (vector_set,
                                    [&v] (const V& e) { return v.partial_order (e).leq (); });
      }

      [[nodiscard]] auto size () const { return vector_set.size (); }

      bool insert (V&& v) {
        bool must_remove = false;

        // This is like remove_if, but allows breaking.
        auto result = vector_set.begin ();
        auto end = vector_set.end ();

        for (auto it = result; it != end; ++it) {
          auto res = v.partial_order (*it);
          if (not must_remove and res.leq ()) {  // v is dominated.
            // if must_remove is true, since we started with an antichain,
            // it's not possible that res.leq () holds.  Hence we don't check for
            // leq if must_remove is true.
            return false;
          }
          if (res.geq ()) {     // v dominates *it
            must_remove = true; /* *it should be removed */
          }
          else {               // *it needs to be kept
            if (result != it)  // This can be false only on the first element.
              *result = std::move (*it);
            ++result;
          }
        }

        if (result != vector_set.end ())
          vector_set.erase (result, vector_set.end ());
        vector_set.push_back (std::move (v));
        return true;
      }

      void union_with (vector_backed&& other) {
        for (auto&& e : other.vector_set)
          insert (std::move (e));
      }

      void intersect_with (const vector_backed& other) {
        vector_backed intersection;
        bool smaller_set = false;

        for (const auto& x : vector_set) {
          bool dominated = false;

          for (auto& y : other.vector_set) {
            V v = x.meet (y);
            if (v == x)
              dominated = true;
            intersection.insert (std::move (v));
            if (dominated)
              break;
          };
          // If x wasn't <= an element in other, then x is not in the
          // intersection, thus the set is updated.
          smaller_set or_eq not dominated;
        }

        if (smaller_set)
          this->vector_set = std::move (intersection.vector_set);
      }

      template <typename F>
      vector_backed apply (const F& lambda) const {
        vector_backed res;
        for (const auto& el : vector_set)
          res.insert (lambda (el));
        return res;
      }

      [[nodiscard]] auto begin () const { return vector_set.begin (); }
      [[nodiscard]] auto end () const { return vector_set.end (); }

      [[nodiscard]] auto& get_backing_vector () { return vector_set; }
      [[nodiscard]] const auto& get_backing_vector () const { return vector_set; }

      template <Vector V2>
      friend class vector_or_kdtree_backed;
  };

  template <Vector V>
  inline std::ostream& operator<< (std::ostream& os, const vector_backed<V>& f) {
    for (auto&& el : f)
      os << el << std::endl;

    return os;
  }
}
