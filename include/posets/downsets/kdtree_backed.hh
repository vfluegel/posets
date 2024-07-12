#pragma once

#include <cassert>
#include <iostream>
#include <cmath>
#include <memory>
#include <set>
#include <vector>

#include <posets/utils/kdtree.hh>

namespace posets::downsets {
  // Forward definition for the operator<<s.
  template <typename>
  class kdtree_backed;

  template <typename Vector>
  std::ostream& operator<<(std::ostream& os, const kdtree_backed<Vector>& f);

  // Another forward def to have friend status
  template <typename Vector>
  class vector_or_kdtree_backed;

  // Finally the actual class definition
  template <typename Vector>
  class kdtree_backed {
    private:
      utils::kdtree<Vector> tree;

      template <typename V>
      friend std::ostream& operator<<(std::ostream& os, const kdtree_backed<V>& f);

      struct proj {
          const Vector& operator() (const Vector* pv) { return *pv; }
          Vector&& operator() (Vector*&& pv)          { return std::move (*pv); }
      };

      void inline resetTree (std::vector<Vector>&& elements) noexcept {
        // we first move everything into a set to remove dupes
        std::set<Vector> elemset;
        for (auto it = elements.begin (); it != elements.end ();)
          elemset.insert (std::move (elements.extract (it++).value ()));
        // and then back into a vector
        std::vector<Vector> noreps;
        noreps.reserve (elemset.size ());
        for (auto it = elemset.begin (); it != elemset.end ();)
          noreps.push_back (std::move (elemset.extract (it++).value ()));

        // now, we can make a tree out of the set to eliminate dominated
        // elements
        this->tree = utils::kdtree<Vector> (std::move (noreps));

        auto antichain = std::vector<Vector*> ();
        antichain.reserve (noreps.size ());
        for (Vector& e : this->tree)
          if (!this->tree.dominates (e, true))
            antichain.push_back (&e);

        this->tree = utils::kdtree<Vector> (
          std::move (antichain), proj ());
        assert (this->tree.is_antichain ());
      }

    public:
      typedef Vector value_type;

      kdtree_backed () = delete;

      kdtree_backed (std::vector<Vector>&& elements) noexcept {
        resetTree (elements);
      }

      kdtree_backed (Vector&& e) :
        tree ( std::array<Vector, 1> { std::move (e) } )
      {}

      template <typename F>
      auto apply (const F& lambda) const {
        const auto& backing_vector = this->tree.vector_set;
        std::vector<Vector> ss;
        ss.reserve (backing_vector.size ());

        for (const auto& v : backing_vector)
          ss.push_back (lambda (v));

        return kdtree_backed (std::move (ss));
      }

      kdtree_backed (const kdtree_backed&) = delete;
      kdtree_backed (kdtree_backed&&) = default;
      kdtree_backed& operator= (const kdtree_backed&) = delete;
      kdtree_backed& operator= (kdtree_backed&&) = default;

      bool contains (const Vector& v) const {
        return this->tree.dominates(v);
      }

      // Union in place
      void union_with (kdtree_backed&& other) {
        assert (other.size () > 0);
        std::vector<Vector*> result;
        result.reserve (this->size () + other.size ());
        // for all elements in this tree, if they are not strictly
        // dominated by the other tree, we keep them
        for (auto& e : tree)
          if (!other.tree.dominates (e, true))
            result.push_back (&e);

        // for all elements in the other tree, if they are not dominated
        // (not necessarily strict) by this tree, we keep them
        for (auto& e : other.tree)
          if (not this->tree.dominates (e))
            result.push_back (&e);
        assert (result.size () > 0);
        this->tree = utils::kdtree<Vector> (
          std::move (result), proj ());
        assert (this->tree.is_antichain ());
      }

      // Intersection in place
      void intersect_with (const kdtree_backed& other) {
        std::vector<Vector> intersection;
        bool smaller_set = false;

        for (auto& x : tree) {
          assert (x.size () > 0);

          // If x is part of the set of all meets, then x will dominate the
          // whole list! So we use this to short-circuit the computation: we
          // first check whether x will be there (which happens only if it is
          // itself dominated by some element in other)
          bool dominated = other.tree.dominates (x);
          if (dominated)
            intersection.push_back (x.copy ());
          else
            for (auto& y : other)
              intersection.push_back (x.meet (y));

          // If x wasn't in the set of meets, dominated is false and
          // the set of minima is different than what is in this->tree
          smaller_set |= not dominated;
        }

        // We can skip building trees and all if this->tree is the antichain
        // of minimal elements
        if (!smaller_set)
          return;

        // Worst-case scenario: we do need to build trees
        resetTree (intersection);
      }

      auto size () const {
        return this->tree.size ();
      }

      auto        begin ()       { return this->tree.begin (); }
      const auto  begin () const { return this->tree.begin (); }
      auto        end ()         { return this->tree.end (); }
      const auto  end () const   { return this->tree.end (); }

      friend class vector_or_kdtree_backed<Vector>;
  };

  template <typename Vector>
  inline std::ostream& operator<<(std::ostream& os, const kdtree_backed<Vector>& f) {
    os << f.tree << std::endl;
    return os;
  }
}
