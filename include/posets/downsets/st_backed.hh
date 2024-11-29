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

  public:
    typedef V value_type;

    static std::unique_ptr<utils::sforest<V>> forest;
    static utils::sforest<V>* get_sforest() {
      return st_backed::forest.get();
    }

    st_backed() = delete;

    st_backed (std::vector<V>&& elements, unsigned k=0) 
    {
      if(!st_backed::forest) {
        st_backed::forest = std::make_unique<utils::sforest<V>>(k, elements.begin()->size());
      }
      this->root = get_sforest()->add_vectors(std::move (elements));
    }

    st_backed (V&& v, unsigned k=0) 
    {
      if(!st_backed::forest) {
        st_backed::forest = std::make_unique<utils::sforest<V>>(k, v.size());
      }
      this->root = get_sforest()->add_vectors(std::array<V, 1> { std::move (v) });
    }

    bool contains (const V& v) const {
        return get_sforest()->covers_vector(this->root, v);
    }

    // Union in place
    void union_with (st_backed&& other) {
        size_t newRoot = get_sforest()->st_union(this->root, other.root);
        this->root = newRoot;
    }

    // Intersection in place
    void intersect_with (const st_backed& other) {
        size_t newRoot = get_sforest()->st_intersect(this->root, other.root);
        this->root = newRoot;
    }

    template <typename F>
    auto apply (const F& lambda) const {
      const auto& all_vectors = get_sforest()->get_all(this->root);
      std::vector<V> ss;
      ss.reserve (all_vectors.size ());

      for (const auto& v : all_vectors)
        ss.push_back (lambda (v));

      return st_backed (std::move (ss));
    }
  };

template <Vector V>
std::unique_ptr<utils::sforest<V>> st_backed<V>::forest = nullptr;

  template <Vector V>
  inline std::ostream& operator<<(std::ostream& os, const st_backed<V>& f) {
    f.get_sforest()->print_children(f.root, 0);
    os << std::endl;
    return os;
  }
}
