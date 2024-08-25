#pragma once

#include <iostream>
#include <functional>
#include <map>
#include <type_traits>
#include <cxxabi.h>
#include <string>

#include <posets/concepts.hh>

template <typename SetType>
class test_t;

template <typename ResType>
struct generic_test {
    virtual ResType operator () () = 0;
    virtual ~generic_test () = default;
};

template <typename ResType>
struct test_list {
    static inline std::map<std::string, std::function<ResType ()> > list;
};

template <typename... Types>
struct type_list;

template <template <typename Arg> typename... Types>
struct template_type_list;

static std::set<std::string> set_names, vector_names;

#define typestring(T)                                                   \
  ([] () -> std::string {                                               \
    int _;                                                              \
    char* s = abi::__cxa_demangle (typeid(T).name (), 0, 0, &_);        \
    std::string ret (s);                                                \
    free (s);                                                           \
    return ret;                                                         \
    }) ()

namespace posets::downsets {
  template <typename VecType>
  struct all {};
}
namespace posets::vectors {
  struct all {};
}

// One set, one vector
template <bool CreateAll = true,
          template <typename V> typename SetType, typename VecType>
void register_maker (type_list<VecType>*, SetType<VecType>* = 0) {
  vector_names.insert (typestring (VecType));

  using res_t = decltype(std::declval<test_t<SetType<VecType>>> () ());
  auto ts = typestring (SetType<VecType>);
  auto test = [ts] () {
    std::cout << "[--] running tests for " << ts << std::endl << std::flush;
    if constexpr (std::is_same_v<res_t, void>) {
      test_t<SetType<VecType>> () ();
      std::cout << "\r[OK]" << std::endl;
    }
    else {
      auto res = test_t<SetType<VecType>> () ();
      std::cout << "\r[OK]" << std::endl;
      return res;
    }
  };

  auto& tests = test_list<std::invoke_result_t<decltype(test)>>::list;
  tests[ts] = test;

  if constexpr (CreateAll) {
    auto setallvec = typestring (SetType<VecType>);
    setallvec = setallvec.substr (0, setallvec.find_first_of ('<')) + "<posets::vectors::all>";

    for (auto&& ts : {
        typestring (posets::downsets::all<VecType>),
        setallvec,
        typestring (posets::downsets::all<posets::vectors::all>)
      }) {
      auto prev = tests[ts];
      tests[ts] = [prev, test] () {
        if (prev) prev ();
        return test ();
      };
    }
  }
}

// One set, multiple vectors
template <bool CreateAll = true,
          template <typename V> typename SetType, typename VecType, typename VecType2, typename... VecTypes>
void register_maker (type_list<VecType, VecType2, VecTypes...>*, SetType<VecType>* = 0) {
  auto s = typestring (SetType<VecType>);
  set_names.insert (s.substr (0, s.find_first_of ('<')));
  register_maker<CreateAll, SetType> ((type_list<VecType>*) 0);
  register_maker<CreateAll, SetType> ((type_list<VecType2, VecTypes...>*) 0);
}

// Multiple sets, multiple vectors
template <bool CreateAll = true,
          template <typename V> typename SetType, template <typename V> typename... SetTypes, typename... VecTypes>
void register_maker (type_list<VecTypes...>* v, template_type_list<SetType, SetTypes...>* = 0) {
  register_maker<CreateAll, SetType> (v);
  register_maker<CreateAll, SetTypes...> (v);
}
