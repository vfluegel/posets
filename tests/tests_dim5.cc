#include <cassert>
#include <span>
#include <memory>
#include <ostream>
#include <set>
#include <vector>
#include <string>
#include <type_traits>
#include <cxxabi.h>

#include "test_maker.hh"

#include <posets/concepts.hh>
#include <posets/downsets.hh>
#include <posets/vectors.hh>

size_t posets::vectors::bool_threshold = 128;
size_t posets::vectors::bitset_threshold = 128;

#define il std::initializer_list<char>

template<typename SetType>
struct test_t : public generic_test<void> {
    using VType = typename SetType::value_type;

    std::vector<VType> vvtovv (const std::vector<std::vector<char>>& vv) {
      std::vector<VType> out;
      //out.reserve (vv.size ());

      for (size_t i = 0; i < vv.size (); ++i)
        out.emplace_back(VType (std::move (vv[i]))); //out[i] = VType<VType> (vv[i]);
      return out;
    }

    SetType vec_to_set (std::vector<VType>&& v) {
      return SetType (std::move (v));
    }

    void operator() () {
      
      // Other tests
      std::cout << "Final tests" << std::endl;
      {
        auto tree = vec_to_set (vvtovv ({
              {7, 0, 9, 9, 7},
              {8, 0, 7, 7, 8},
              {8, 0, 9, 9, 8},
              {9, 0, 7, 7, 9}
            }));
        //assert (tree.size () == 4);
        assert (tree.contains (VType (il {0, 0, 0, 0, 0})));
        assert (tree.contains (VType (il {6, 0, 9, 9, 7})));
        assert (tree.contains (VType (il {7, 0, 9, 9, 7})));
      }


      // Full set is so slow, this will never finish.
      if constexpr (std::is_same<SetType, posets::downsets::full_set<VType>>::value)
        return;


      auto F1i = vec_to_set (vvtovv ({
            {7, 0, 9, 9, 7},
            {8, 0, 9, 9, 6},
            {9, 0, 7, 7, 9}
          }));
      auto F = vec_to_set (vvtovv ({
            {7, 0, 9, 9, 7},
            {8, 0, 9, 9, 5},
            {9, 0, 7, 7, 9}
          }));

      std::cout << "Preparing to intersect" << std::endl;
      F.intersect_with (std::move (F1i));

      assert (F.contains (VType (il {8, 0, 9 ,9, 4})));
      assert (F.contains (VType (il {8, 0, 9 ,9, 5})));
      assert (not F.contains (VType (il {8, 0, 9 ,9, 6})));
      assert (F.contains (VType (il {7, 0, 9 ,9, 7})));
      assert (not F.contains (VType (il {7, 0, 9 ,9, 8})));
      assert (F.contains (VType (il {9, 0, 7 ,7, 9})));
      assert (not F.contains (VType (il {9, 0, 7 ,7, 10})));
      
    }

};



void usage (const char* progname) {
  std::cout << "usage: " << progname << " SETTYPE VECTYPE" << std::endl;

  std::cout << "  SETTYPE:\n  - all" << std::endl;
  for (auto&& s : set_names) {
    auto start = s.find_last_of (':') + 1;
    std::cout << "  - " << s.substr (start, s.find_first_of ('<') - start) << std::endl;
  }

  std::cout << "  VECTYPE:\n  - all" << std::endl;
  for (auto&& s : vector_names) {
    auto start = s.find_last_of (':') + 1;
    std::cout << "  - " << s.substr (start, s.find_first_of ('>') - start + 1) << std::endl;
  }

  exit (0);
}

namespace posets::vectors {
  template <typename T>
  using array_backed_fixed = posets::vectors::array_backed<T, 10>;

  template <typename T>
  using array_backed_sum_fixed = posets::vectors::array_backed_sum<T, 10>;

  template <typename T>
  using simd_array_backed_fixed = posets::vectors::simd_array_backed<T, 10>;

  template <typename T>
  using simd_array_ptr_backed_fixed = posets::vectors::simd_array_ptr_backed<T, 10>;

  template <typename T>
  using simd_array_backed_sum_fixed = posets::vectors::simd_array_backed_sum<T, 10>;

  template <typename T>
  using simd_array_ptr_backed_sum_fixed = posets::vectors::simd_array_ptr_backed_sum<T, 10>;

  template <typename T>
  using simd_vector_and_bitset_backed = posets::vectors::X_and_bitset<posets::vectors::simd_vector_backed<T>, 1>;
}

#define DEFINE_VECTOR_NAME(V) template <> struct vector_name<V> { static constexpr auto str = #V; };
#define DEFINE_VECTOR_NAMES(...)
#define VECTOR_TYPES(...) FOR_EACH(DEFINE_VECTOR_NAME, __VA_ARGS__)     \
  using vector_types = type_list<__VA_ARGS__>;

VECTOR_TYPES (posets::vectors::vector_backed<char>,
              posets::vectors::array_backed_fixed<char>,
              posets::vectors::array_backed_sum_fixed<char>,
              posets::vectors::simd_vector_backed<char>,
              posets::vectors::simd_array_backed_fixed<char>,
              posets::vectors::simd_array_ptr_backed_fixed<char>,
              posets::vectors::simd_array_backed_sum_fixed<char>,
              posets::vectors::simd_array_ptr_backed_sum_fixed<char>,
              posets::vectors::simd_vector_and_bitset_backed<char>);

using set_types = template_type_list<//posets::downsets::full_set, ; too slow.
  posets::downsets::sharingtree_backed,
  posets::downsets::kdtree_backed,
  posets::downsets::vector_or_kdtree_backed,
  posets::downsets::set_backed,
  posets::downsets::vector_backed,
  posets::downsets::vector_backed_bin,
  posets::downsets::vector_backed_one_dim_split,
  posets::downsets::vector_backed_one_dim_split_intersection_only>;



int main(int argc, char* argv[]) {
  register_maker ((vector_types*) 0, (set_types*) 0);

  if (argc != 3)
    usage (argv[0]);

  auto implem = std::string ("posets::downsets::") + argv[1]
    + "<posets::vectors::" + argv[2] + ">";

  try {
    posets::vectors::bool_threshold = 128;
    posets::vectors::bitset_threshold = 128;
    auto& tests = test_list<void>::list;
    tests[implem] ();
    return 0;
  } catch (std::bad_function_call& e) {
    std::cout << "error: no such implem: " << implem << std::endl;
    return 1;
  }
}
