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

    assert(f.cover_vector(idcs, std::move(vvtovv(data))));

    return 0;
}
