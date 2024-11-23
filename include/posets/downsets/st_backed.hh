#pragma once

#include <cassert>
#include <iostream>
#include <cmath>
#include <memory>
#include <set>
#include <vector>

#include <posets/concepts.hh>
#include <posets/utils/sforest.hh>

namespace posets::downsets {
  // Forward definition for the operator<<s.
  template <Vector>
  class st_backed;

  template <Vector V>
  std::ostream& operator<<(std::ostream& os, const st_backed<V>& f);

  template <Vector V> class st_backed {
  private:
    template <Vector V2>
    friend std::ostream& operator<<(std::ostream& os, const st_backed<V2>& f);

    int k;
    size_t dim;

    size_t root{};
    static utils::sforest<V> get_sforest() {
        static utils::sforest<V> f(k, dim);
        return f;
    }

  public:
    st_backed() = delete;

    st_backed (std::vector<V>&& elements) 
        : k {}
        , dim {elements.begin().size()} 
    {
        get_sforest().add_vectors(std::move (elements));
    }

    st_backed (V&& v) 
        : k {}
        , dim {v.size()} 
    {
        get_sforest().add_vectors(std::array<V, 1> { std::move (elements) });
    }

    bool contains (const V& v) const {
        return this->f.cover_vector(this->root, v);
    }

    // Union in place
    void union_with (st_backed&& other) {
        size_t newRoot = f.st_union(this->root, other.root);
        this->root = newRoot;
    }

    // Intersection in place
    void intersect_with (const kdtree_backed& other) {
        size_t newRoot = f.st_intersect(this->root, other.root);
        this->root = newRoot;
    }
  };

  template <Vector V>
  inline std::ostream& operator<<(std::ostream& os, const st_backed<V>& f) {
    f.get_sforest().print_children(f.root, 0);
    os << std::endl;
    return os;
  }
}