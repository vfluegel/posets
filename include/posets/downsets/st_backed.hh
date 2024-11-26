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
    static std::unique_ptr<utils::sforest<V>> forest{};

  public:
    static utils::sforest<V>* get_sforest() {
      return forest.get();
    }

    st_backed() = delete;

    st_backed (std::vector<V>&& elements, unsigned k=0) 
    {
      if(!forest) {
        forest = std::make_unique<utils::sforest<V>>(k, elements.begin().size());
      }
      get_sforest()->add_vectors(std::move (elements));
    }

    st_backed (V&& v, unsigned k=0) 
    {
      if(!forest) {
        forest = std::make_unique<utils::sforest<V>>(k, v.size());
      }
      get_sforest()->add_vectors(std::array<V, 1> { std::move (v) });
    }

    bool contains (const V& v) const {
        return get_sforest()->cover_vector(this->root, v);
    }

    // Union in place
    void union_with (st_backed&& other) {
        size_t newRoot = get_sforest()->st_union(this->root, other.root);
        this->root = newRoot;
    }

    // Intersection in place
    void intersect_with (const kdtree_backed& other) {
        size_t newRoot = get_sforest()->st_intersect(this->root, other.root);
        this->root = newRoot;
    }
  };

  template <Vector V>
  inline std::ostream& operator<<(std::ostream& os, const st_backed<V>& f) {
    f.get_sforest()->print_children(f.root, 0);
    os << std::endl;
    return os;
  }
}
