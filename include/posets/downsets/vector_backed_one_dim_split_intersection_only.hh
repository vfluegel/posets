#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <vector>

#include <posets/concepts.hh>
#include <posets/utils/vector_mm.hh>

namespace posets::downsets {
  template <Vector V>
  class vector_backed_one_dim_split_intersection_only {
      using self = vector_backed_one_dim_split_intersection_only;

    public:
      using value_type = V;

      vector_backed_one_dim_split_intersection_only (V&& v) { insert (std::move (v)); }

      vector_backed_one_dim_split_intersection_only (std::vector<V>&& elements) noexcept {
        assert (not elements.empty ());
        for (auto&& e : elements)
          insert (std::move (e));
      }

    private:
      vector_backed_one_dim_split_intersection_only () = default;

    public:
      vector_backed_one_dim_split_intersection_only (const self&) = delete;
      vector_backed_one_dim_split_intersection_only (self&&) = default;
      self& operator= (self&&) = default;
      self& operator= (const self&) = delete;

      bool operator== (const self& other) = delete;

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

      void union_with (self&& other) {
        for (auto&& e : other.vector_set)
          insert (std::move (e));
      }

      template <typename RefVec>
      class disregard_first_component {
        public:
          bool operator() (const RefVec& v1, const RefVec& v2) const {
            utils::vector_mm<typename V::value_type> v (v1.get ().size ());
            v.reserve (v.size ());

            v2.get ().to_vector (v);
            v[0] = v1.get ()[0];
            v.resize (v2.get ().size ());
            return v1.get () < V (v);
          }
      };

      void intersect_with (const self& other) {
        self intersection;
        bool smaller_set = false;

        // split_cache maps first-dim p to:
        //   - the set of vectors that have a bigger or equal to p
        //   - the set of other vectors
        // then we go through all vectors in vector_set:
        //    for all vectors in split_cache[x[0]].first, compute the meets, and
        //              put them intersection; if x is dominated, leave.
        //    if x is not dominated, do the same for split_cache[x[0]].second

        // Different wording: we compare x in vector_set with all of
        // other.vector_set, prioritazing the vectors in other.vector_set that
        // have a larger first component (hence can potentially dominate x).

        using cache_red_dim = std::set<std::reference_wrapper<const V>,
                                       disregard_first_component<std::reference_wrapper<const V>>>;
        using vector_of_vectors = std::vector<std::reference_wrapper<const V>>;
        std::map<int, std::pair<cache_red_dim, vector_of_vectors>> split_cache;

        for (const auto& x : vector_set) {
          bool dominated = false;

          auto& cv = ([&x, &split_cache, &other] () -> auto& {
            try {
              return split_cache.at (x[0]);
            } catch (...) {
              auto& cv = split_cache[x[0]];
              for (const auto& y : other.vector_set)
                if (y[0] >= x[0])
                  cv.first.insert (std::ref (y));
                else
                  cv.second.push_back (std::ref (y));
              return cv;
            }
          }) ();

          auto meet = [&] (std::reference_wrapper<const V> y) {
            V&& v = x.meet (y.get ());
            if (v == x)
              dominated = true;
            intersection.insert (std::move (v));
            return not dominated;
          };

          static_cast<void> (std::all_of (cv.first.begin (), cv.first.end (), meet));
          if (not dominated)
            static_cast<void> (std::all_of (cv.second.begin (), cv.second.end (), meet));

          // If x wasn't <= an element in other, then x is not in the
          // intersection, thus the set is updated.
          smaller_set or_eq not dominated;
        }

        if (smaller_set)
          this->vector_set = std::move (intersection.vector_set);
      }

      template <typename F>
      self apply (const F& lambda) const {
        self res;
        for (const auto& el : vector_set)
          res.insert (lambda (el));
        return res;
      }

      [[nodiscard]] auto& get_backing_vector () { return vector_set; }

      [[nodiscard]] const auto& get_backing_vector () const { return vector_set; }

      [[nodiscard]] auto begin () { return vector_set.begin (); }
      [[nodiscard]] auto begin () const { return vector_set.begin (); }
      [[nodiscard]] auto end () { return vector_set.end (); }
      [[nodiscard]] auto end () const { return vector_set.end (); }

    private:
      std::vector<V> vector_set;
  };

  template <Vector V>
  inline std::ostream& operator<< (std::ostream& os,
                                   const vector_backed_one_dim_split_intersection_only<V>& f) {
    for (auto&& el : f)
      os << el << std::endl;

    return os;
  }
}
