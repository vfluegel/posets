#pragma once

#include <posets/concepts.hh>
#include <posets/downsets.hh>

#include <posets/vectors/vector_backed.hh>
#include <posets/vectors/array_backed.hh>
#include <posets/vectors/array_ptr_backed.hh>
#include <posets/vectors/array_backed_sum.hh>
#include <posets/vectors/simd_vector_backed.hh>
#include <posets/vectors/simd_array_backed.hh>
#include <posets/vectors/X_and_bitset.hh>

namespace std {
  template <posets::Vector V>
  inline
  std::ostream& operator<<(std::ostream& os, const V& v) {  return v.print (os); }
}
