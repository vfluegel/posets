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
  
  template <Vector V> class sharingtree_backed {
  private:
    size_t dim;

    size_t root{};

  public:
    typedef V value_type;
    
    std::vector<V> vector_set;
    static std::unique_ptr<utils::sharingforest<V>> forest;
    static utils::sharingforest<V>* sharingforest() {
      return sharingtree_backed::forest.get();
    }

    sharingtree_backed() = delete;
    sharingtree_backed (const sharingtree_backed&) = delete;
    sharingtree_backed (sharingtree_backed&&) = default;
    sharingtree_backed& operator= (const sharingtree_backed&) = delete;
    sharingtree_backed& operator= (sharingtree_backed&&) = default;

    sharingtree_backed (std::vector<V>&& elements) 
    {
      if(!sharingtree_backed::forest) {
        sharingtree_backed::forest = std::make_unique<utils::sharingforest<V>>(elements.begin()->size());
      }
      this->root = sharingforest()->add_vectors(std::move (elements));
      this->vector_set = sharingforest()->get_all(this->root);
    }

    sharingtree_backed (V&& v) 
    {
      if(!sharingtree_backed::forest) {
        sharingtree_backed::forest = std::make_unique<utils::sharingforest<V>>(v.size());
      }
      this->root = sharingforest()->add_vectors(std::array<V, 1> { std::move (v) });
      this->vector_set = sharingforest()->get_all(this->root);
    }

    auto size () const {
      return this->vector_set.size();
    }
    auto        begin ()       { return this->vector_set.begin (); }
    const auto  begin () const { return this->vector_set.begin (); }
    auto        end ()         { return this->vector_set.end (); }
    const auto  end () const   { return this->vector_set.end (); }

    bool contains (const V& v) const {
        return sharingforest()->covers_vector(this->root, v);
    }

    // Union in place
    void union_with (sharingtree_backed&& other) {
        size_t newRoot = sharingforest()->st_union(this->root, other.root);
        this->root = newRoot;
        this->vector_set = sharingforest()->get_all(this->root);
    }

    // Intersection in place
    void intersect_with (const sharingtree_backed& other) {
        size_t newRoot = sharingforest()->st_intersect(this->root, other.root);
        this->root = newRoot;
        this->vector_set = sharingforest()->get_all(this->root);
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
std::unique_ptr<utils::sharingforest<V>> sharingtree_backed<V>::forest = nullptr;

  template <Vector V>
  inline std::ostream& operator<<(std::ostream& os, const sharingtree_backed<V>& f) {
    for (auto &&el : f.vector_set)
      os << el << std::endl;
    return os;
  }
}
