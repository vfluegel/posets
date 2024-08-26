#pragma once

namespace posets::vectors {
  /// Conditional member when has_sum is set in @a simd_array_backed_.
  template <bool has_sum>
  struct sum_member { int sum; };

  template<>
  struct sum_member<false> {};

  /// Conditional member when embeds_data is unset in @a simd_array_backed_.
  template <bool embeds_data, typename Data>
  struct malloc_member {
      static inline boost::object_pool<Data> malloc {8192};
  };

  template <typename Data>
  struct malloc_member<true, Data> {};
}
