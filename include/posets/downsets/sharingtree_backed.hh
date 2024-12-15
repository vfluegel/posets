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
    std::shared_ptr<utils::sharingforest<V>> forest;
    std::vector<V> vector_set;
    static std::map<size_t, std::shared_ptr<utils::sharingforest<V>>> forest_map;

    void init_forest(size_t dimkey) {
      auto res = sharingtree_backed::forest_map.find (dimkey);
      if (res != sharingtree_backed::forest_map.end ()) {
        this->forest = res->second;
      } else {
        this->forest = std::make_shared<utils::sharingforest<V>> (dimkey);
        sharingtree_backed::forest_map.emplace(dimkey, this->forest);
      }
    }

  public:
    typedef V value_type;   

    sharingtree_backed() = delete;
    sharingtree_backed (const sharingtree_backed&) = delete;
    sharingtree_backed (sharingtree_backed&&) = default;
    sharingtree_backed& operator= (const sharingtree_backed&) = delete;
    sharingtree_backed& operator= (sharingtree_backed&&) = default;

    sharingtree_backed (std::vector<V>&& elements) 
    {
      init_forest (elements.begin()->size());
      this->root = this->forest->add_vectors(std::move (elements));
      this->vector_set = this->forest->get_all(this->root);
    }

    sharingtree_backed (V&& v) 
    {
      init_forest (v.size ());
      this->root = this->forest->add_vectors(std::array<V, 1> { std::move (v) });
      this->vector_set = this->forest->get_all(this->root);
    }

    auto size () const {
      return this->vector_set.size();
    }
    auto        begin ()       { return this->vector_set.begin (); }
    const auto  begin () const { return this->vector_set.begin (); }
    auto        end ()         { return this->vector_set.end (); }
    const auto  end () const   { return this->vector_set.end (); }

    bool contains (const V& v) const {
        return this->forest->covers_vector(this->root, v);
    }

    // Union in place
    void union_with (sharingtree_backed&& other) {
        size_t newRoot = this->forest->st_union(this->root, other.root);
        this->root = newRoot;
        this->vector_set = this->forest->get_all(this->root);
    }

    // Intersection in place
    void intersect_with (const sharingtree_backed& other) {
        size_t newRoot = this->forest->st_intersect(this->root, other.root);
        this->root = newRoot;
        this->vector_set = this->forest->get_all(this->root);
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
  std::map<size_t, std::shared_ptr<utils::sharingforest<V>>> sharingtree_backed<V>::forest_map;

  template <Vector V>
  inline std::ostream& operator<<(std::ostream& os, const sharingtree_backed<V>& f) {
    for (auto &&el : f.vector_set)
      os << el << std::endl;
    return os;
  }
}
