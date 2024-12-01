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

#include <posets/downsets.hh>
#include <posets/vectors.hh>

size_t posets::vectors::bool_threshold = 128;
size_t posets::vectors::bitset_threshold = 128;

template<class T, class = void>
struct has_insert : std::false_type {};

template <class T>
struct has_insert<T, std::void_t<decltype(std::declval<T>().insert (std::declval<typename T::value_type> ()))>> : std::true_type {};

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
      if constexpr (has_insert<SetType>::value) {
        SetType set (std::move (v[0]));
        for (size_t i = 1; i < v.size (); ++i)
          set.insert (std::move (v[i]));
        return set;
      }
      else
        return SetType (std::move (v), 11);
    }

    void operator() () {
      
      auto F = vec_to_set (vvtovv ({
            {0, 7, 0, 0, 9, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 8, 0, 0, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-1, 8, -1, 0, 9, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-1, 8, -1, 0, -1, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-1, 7, -1, 0, -1, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-1, 8, -1, 0, -1, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-1, 8, -1, 0, -1, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {-1, 9, -1, 0, -1, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
          }));
      assert (F.contains (VType (il {-1, 9, -1, 0, -1, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})));
    
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
  posets::downsets::st_backed,
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
