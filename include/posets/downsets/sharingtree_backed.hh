#pragma once

#include <cassert>
#include <iostream>
#include <cmath>
#include <memory>
#include <set>
#include <vector>

#include <posets/concepts.hh>
#include <posets/utils/sharingforest.hh>

namespace posets::downsets {
  // Forward definition for the operator<<s.
  template <Vector>
  class sharingtree_backed;

  template <Vector V>
  std::ostream& operator<<(std::ostream& os, const sharingtree_backed<V>& f);

  template <Vector V> class sharingtree_backed {
  private:
    template <Vector V2>
    friend std::ostream& operator<<(std::ostream& os, const sharingtree_backed<V2>& f);

    int k;
    size_t dim;

    size_t root{};

  public:
    typedef V value_type;

    static std::unique_ptr<utils::sharingforest<V>> forest;
    static utils::sharingforest<V>* sharingforest() {
      return sharingtree_backed::forest.get();
    }

    sharingtree_backed() = delete;

    sharingtree_backed (std::vector<V>&& elements, unsigned k=0) 
    {
      if(!sharingtree_backed::forest) {
        sharingtree_backed::forest = std::make_unique<utils::sharingforest<V>>(k, elements.begin()->size());
      }
      this->root = sharingforest()->add_vectors(std::move (elements));
    }

    sharingtree_backed (V&& v, unsigned k=0) 
    {
      if(!sharingtree_backed::forest) {
        sharingtree_backed::forest = std::make_unique<utils::sharingforest<V>>(k, v.size());
      }
      this->root = sharingforest()->add_vectors(std::array<V, 1> { std::move (v) });
    }

    bool contains (const V& v) const {
        return sharingforest()->covers_vector(this->root, v);
    }

    // Union in place
    void union_with (sharingtree_backed&& other) {
        size_t newRoot = sharingforest()->st_union(this->root, other.root);
        this->root = newRoot;
    }

    // Intersection in place
    void intersect_with (const sharingtree_backed& other) {
        size_t newRoot = sharingforest()->st_intersect(this->root, other.root);
        this->root = newRoot;
    }

    template <typename F>
    auto apply (const F& lambda) const {
      const auto& all_vectors = sharingforest()->get_all(this->root);
      std::vector<V> ss;
      ss.reserve (all_vectors.size ());

      for (const auto& v : all_vectors)
        ss.push_back (lambda (v));

      return sharingtree_backed (std::move (ss));
    }
  };

template <Vector V>
std::unique_ptr<utils::sharingforest<V>> sharingtree_backed<V>::forest = nullptr;

  template <Vector V>
  inline std::ostream& operator<<(std::ostream& os, const sharingtree_backed<V>& f) {
    f.sharingforest()->print_children(f.root, 0);
    os << std::endl;
    return os;
  }
}
