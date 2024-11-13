#include <vector>

#include <posets/utils/sforest.hh>
#include <posets/vectors.hh>

namespace utils = posets::utils;

using VType = posets::vectors::vector_backed<char>;

std::vector<VType> vvtovv (const std::vector<std::vector<char>>& vv) {
  std::vector<VType> out;
  for (size_t i = 0; i < vv.size (); ++i)
    out.emplace_back(VType (std::move (vv[i])));
  return out;
}

int main(int argc, char const *argv[])
{
    std::vector<std::vector<char>> data{{6, 3, 2}, {5, 5, 4}, {2, 6, 2}};

    utils::sforest<VType> f{10, 3};
    auto idcs = f.add_vectors(std::move(vvtovv(data)));

    // Test adding a second tree to the forest
    std::vector<std::vector<char>> data2{{7, 5, 3}, {4, 8, 4}, {2, 5, 6}};
    auto idcs2 = f.add_vectors(std::move(vvtovv(data2)));

    std::vector<char> v = {2, 2, 2};
    assert(f.cover_vector(idcs, VType (std::move(v))));

    return 0;
}
