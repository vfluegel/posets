#pragma once

#include <cstring>
#include <cassert>
#include <experimental/simd>
#include <iostream>

#include <posets/concepts.hh>
#include <posets/vectors/simd_po_res.hh>
#include <posets/utils/simd_traits.hh>

#include <posets/vectors/sum_ptr_helpers.hh>

namespace posets::vectors {
  template <typename T, size_t nsimds, bool has_sum, bool embeds_data>
  class simd_array_backed_;

  template <typename T, size_t K>
  using simd_array_backed = simd_array_backed_<T, utils::simd_traits<T>::nsimds (K), false, true>;

  template <typename T, size_t K>
  using simd_array_backed_sum = simd_array_backed_<T, utils::simd_traits<T>::nsimds (K), true, true>;

  template <typename T, size_t K>
  using simd_array_ptr_backed = simd_array_backed_<T, utils::simd_traits<T>::nsimds (K), false, false>;

  template <typename T, size_t K>
  using simd_array_ptr_backed_sum = simd_array_backed_<T, utils::simd_traits<T>::nsimds (K), true, false>;

  template <typename T, size_t nsimds, bool has_sum, bool embeds_data>
  class simd_array_backed_ :
      private sum_member<has_sum>,
      private malloc_member<embeds_data,
                            std::array<typename utils::simd_traits<T>::fssimd, nsimds>> {

      using self = simd_array_backed_;
      using traits = utils::simd_traits<T>;
      static const auto simd_size = traits::simd_size;
      using data_t = std::array<typename traits::fssimd, nsimds>;

    private:
      data_t& get_data () {
        if constexpr (embeds_data)
          return data;
        else
          return *data;
      }

      const data_t& get_data () const {
        if constexpr (embeds_data)
          return data;
        else
          return *data;
      }

    public:
      using value_type = T;

      simd_array_backed_ (size_t k) : k {k} {
        if constexpr (not embeds_data)
          data = this->malloc.construct ();
        if constexpr (has_sum)
          this->sum = 0;
        assert (utils::simd_traits<T>::nsimds (k) == nsimds);
        if (k % simd_size != 0)
          get_data ().back () = 0;
      }

      simd_array_backed_ (std::span<const T> v) : simd_array_backed_ (v.size ()) {
        if constexpr (has_sum) {
          this->sum = 0;
          for (auto&& c : v)
            this->sum += c;
        }
        if (k % simd_size != 0)
          get_data ().back () = 0;
        std::memcpy ((char*) get_data ().data (), (char*) v.data (), v.size () * sizeof (T));
      }

      simd_array_backed_ () = delete;
      simd_array_backed_ (const self& other) = delete;
      simd_array_backed_ (self&& other) : k {other.k}, data {other.data} {
        if constexpr (not embeds_data) other.data = nullptr;
        if constexpr (has_sum) this->sum = other.sum;
      }

      ~simd_array_backed_ () {
        if constexpr (not embeds_data) if (data) this->malloc.destroy (data);
      }

      // explicit copy operator
      simd_array_backed_ copy () const {
        if constexpr (embeds_data) {
          auto res = simd_array_backed_ (k);
          res.data = data;
          if constexpr (has_sum) res.sum = this->sum;
          return res;
        }
        else {
          auto res = simd_array_backed_ (std::span ((T*) data, k));
          if constexpr (has_sum) res.sum = this->sum;
          return res;
        }
      }

      self& operator= (self&& other) {
        if constexpr (not embeds_data)
          if (data) this->malloc.destroy (data);
        k = other.k;
        data = other.data;
        if constexpr (not embeds_data) other.data = nullptr;
        if constexpr (has_sum) this->sum = other.sum;

        return *this;
      }

      self& operator= (const self& other) = delete;

      static constexpr size_t capacity_for (size_t elts) {
        return nsimds * simd_size;
      }

      void to_vector (std::span<T> v) const {
        std::memcpy ((char*) v.data (), (char*) get_data ().data (), k * sizeof (T));
      }


      inline auto partial_order (const self& rhs) const {
        return simd_po_res (*this, rhs);
      }

      // Used by Sets, should be a total order.  Do not use.
      bool operator< (const self& rhs) const {
        for (size_t i = 0; i < nsimds; ++i) {
          auto lhs_lt_rhs = get_data ()[i] < rhs.get_data ()[i];
          auto rhs_lt_lhs = rhs.get_data ()[i] < get_data ()[i];
          auto p1 = find_first_set (lhs_lt_rhs);
          auto p2 = find_first_set (rhs_lt_lhs);
          if (p1 == p2)
            continue;
          return (p1 < p2);
        }
        return false;
      }

      bool operator== (const self& rhs) const {
        if constexpr (has_sum)
          if (this->sum != rhs.sum)
            return false;
        // Trust memcmp to DTRT
        return std::memcmp ((char*) rhs.get_data ().data (), (char*) get_data ().data (), k * sizeof (T)) == 0;
      }

      bool operator!= (const self& rhs) const {
        if constexpr (has_sum)
          if (this->sum != rhs.sum)
            return true;
        // Trust memcmp to DTRT
        return std::memcmp ((char*) rhs.get_data ().data (), (char*) get_data ().data (), k * sizeof (T)) != 0;
      }

      self meet (const self& rhs) const {
        auto res = self (k);
        if constexpr (not embeds_data) res.data = this->malloc.construct ();

        for (size_t i = 0; i < nsimds; ++i) {
          res.get_data ()[i] = std::experimental::min (get_data ()[i], rhs.get_data ()[i]);

          // The following should NOT be used since this can lead to overflows over char:
          //   res.sum += std::experimental::reduce (res.get_data ()[i]);
          // instead, we manually loop through:
          if constexpr (has_sum)
            for (size_t j = 0; j < simd_size; ++j)
              res.sum += res.get_data ()[i][j];
        }

        return res;
      }

      auto size () const {
        return k;
      }

      auto& print (std::ostream& os) const
      {
        os << "{ ";
        for (size_t i = 0; i < k; ++i)
          os << (int) (*this)[i] << " ";
        os << "}";
        return os;
      }

      // Should be used sparingly.
      int operator[] (size_t i) const {
        assert (i / simd_size < nsimds);
        return get_data ()[i / simd_size][i % simd_size];
      }

      auto bin () const {
        if constexpr (has_sum)
          return std::abs (this->sum) / k;
        else
          return std::abs ((*this)[0]) / k;
      }

    private:
      friend simd_po_res<self>;
      size_t k;
      std::conditional_t<embeds_data, data_t, data_t*> data;
  };


  template <typename T>
  concept SimdArrayBacked = requires (T&& u) { new simd_array_backed_ (std::move (u)); };

  template <template <typename, auto...> typename S, typename T>
  requires SimdArrayBacked<S<T, 64>>
  struct traits<S, T> {
      static constexpr auto capacity_for (size_t elts) {
        return utils::simd_traits<T>::capacity_for (elts);
      }
  };

  static_assert (Vector<simd_array_backed<int, 128>>);
  static_assert (Vector<simd_array_backed_sum<int, 128>>);
  static_assert (Vector<simd_array_ptr_backed<int, 128>>);
  static_assert (Vector<simd_array_ptr_backed_sum<int, 128>>);
}
