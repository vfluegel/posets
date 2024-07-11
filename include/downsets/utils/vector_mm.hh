#pragma once

// https://stackoverflow.com/questions/8456236/how-is-a-vectors-data-aligned

#include <vector>
#include <boost/align/aligned_allocator.hpp>
#include <downsets/utils/simd_traits.hh>

namespace downsets::utils {
  template <typename T>
  using vector_mm = std::vector<T, boost::alignment::aligned_allocator<T, simd_traits<T>::alignment ()>>;
}