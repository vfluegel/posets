#pragma once

#include <cassert>
#include <iostream>
#include <span>

namespace posets::vectors {
  // What's the multiple of T's we store.  This is used to speed up compilation
  // and reduce program size.
#define T_PER_UNIT 8

  template <typename T, size_t Units>
  class array_backed_;

  template <typename T, size_t K>
  using array_backed = array_backed_<T, (K + T_PER_UNIT - 1) / T_PER_UNIT>;

  template <typename T, size_t Units>
  class array_backed_ : public std::array<T, Units * T_PER_UNIT> {
      using self = array_backed_<T, Units>;
      using base = std::array<T, Units * T_PER_UNIT>;

    public:
      array_backed_ (size_t k) : k {k} { assert (k <= Units * T_PER_UNIT); }

      array_backed_ (std::span<const T> v) : array_backed_ (v.size ()) {
        std::copy (v.begin (), v.end (), this->data ());
        if (Units * T_PER_UNIT > k)
          std::fill (&(*this)[k], this->end (), 0);
      }

      size_t size () const { return k; }

    private:
      array_backed_ (const self& other) = default;

    public:
      array_backed_ (self&& other) = default;

      // explicit copy operator
      array_backed_ copy () const {
        auto res = array_backed_ (*this);
        return res;
      }

      self& operator= (self&& other) {
        assert (other.k == k);
        base::operator= (std::move (other));
        return *this;
      }

      self& operator= (const self& other) = delete;

      static constexpr size_t capacity_for (size_t elts) {
        return elts;
      }

      void to_vector (std::span<T> v) const {
        assert (v.size () >= k);
        std::copy (this->begin (), &(*this)[k], v.data ());
      }

      class po_res {
        public:
          po_res (const self& lhs, const self& rhs) {
            bgeq = true;
            bleq = true;
            for (size_t i = 0; i < lhs.k; ++i) {
              auto diff = lhs[i] - rhs[i];
              if (bgeq)
                bgeq = (diff >= 0);
              if (bleq)
                bleq = (diff <= 0);
              if (not bgeq and not bleq)
                break;
            }
          }

          inline bool geq () {
            return bgeq;
          }

          inline bool leq () {
            return bleq;
          }
        private:
          bool bgeq, bleq;
      };

      inline auto partial_order (const self& rhs) const {
        assert (rhs.k == k);

        return po_res (*this, rhs);
      }

      self meet (const self& rhs) const {
        const auto ub = std::min (k, Units * T_PER_UNIT);
        auto res = self (k);
        assert (rhs.k == k);
        for (size_t i = 0; i < ub; ++i)
          res[i] = std::min ((*this)[i], rhs[i]);

        if (k < Units * T_PER_UNIT)
          std::fill (&res[k], res.end (), 0);

        return res;
      }

      std::ostream& print (std::ostream& os) const {
        os << "{ ";
        for (size_t i = 0; i < this->size (); ++i)
          os << (int) (*this)[i] << " ";
        os << "}";
        return os;
      }

    private:
      const size_t k;
  };

  static_assert (Vector<array_backed<int, 128>>);
}
