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
    auto root1 = f.add_vectors(std::move(vvtovv(data)));

    // Test adding a second tree to the forest
    std::vector<std::vector<char>> data2{{7, 5, 3}, {4, 8, 4}, {2, 5, 6}};
    auto root2 = f.add_vectors(std::move(vvtovv(data2)));

    std::cout << f << std::endl;
    f.print_children(root1, 0);
    std::cout << std::endl;
    std::vector<char> v = {2, 2, 2};
    assert(f.cover_vector(root1, VType (std::move(v))));
    std::vector<char> w = {9, 4, 8};
    assert(!f.cover_vector(root1, VType (std::move(w))));

    auto iRoot = f.st_intersect(root1, root2);
    f.print_children(iRoot, 0);
    std::cout << std::endl;

    return 0;
}
