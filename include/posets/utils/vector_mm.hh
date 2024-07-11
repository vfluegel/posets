#pragma once

// https://stackoverflow.com/questions/8456236/how-is-a-vectors-data-aligned

#include <vector>
#include <boost/align/aligned_allocator.hpp>
#include <posets/utils/simd_traits.hh>

namespace posets::utils {
  template <typename T>
  using vector_mm = std::vector<T, boost::alignment::aligned_allocator<T, simd_traits<T>::alignment ()>>;
}