#include <vector>

#include <posets/utils/sforest.hh>
#include <posets/vectors.hh>

namespace utils = posets::utils;

using VType = posets::vectors::vector_backed<char>;

std::vector<VType> vvtovv(const std::vector<std::vector<char>> &vv) {
  std::vector<VType> out;
  for (size_t i = 0; i < vv.size(); ++i)
    out.emplace_back(VType(std::move(vv[i])));
  return out;
}

int main(int argc, char const *argv[]) {
  utils::sforest<VType> f{10, 3};

  // Add some vectors to create a tree
  std::vector<std::vector<char>> data{{6, 3, 2}, {5, 5, 4}, {2, 6, 2}};
  auto idcs = f.add_vectors(std::move(vvtovv(data)));
  std::cout << f << std::endl;
  std::vector<char> v = {2, 2, 2};
  assert(f.covers_vector(idcs, VType(std::move(v))));
  v = {6, 3, 1};
  assert(f.covers_vector(idcs, VType(std::move(v))));
  v = {8, 3, 1};
  assert(not f.covers_vector(idcs, VType(std::move(v))));
  v = {9, 4, 8};
  assert(not f.covers_vector(idcs, VType(std::move(v))));

  // Test adding a second tree to the forest
  data = {{7, 4, 3}, {4, 8, 4}, {2, 5, 6}};
  auto idcs2 = f.add_vectors(std::move(vvtovv(data)));
  std::cout << f << std::endl;

  // We repeat the tests above to check we did not break the previous tree
  v = {2, 2, 2};
  assert(f.covers_vector(idcs, VType(std::move(v))));
  v = {6, 3, 1};
  assert(f.covers_vector(idcs, VType(std::move(v))));
  v = {8, 3, 1};
  assert(not f.covers_vector(idcs, VType(std::move(v))));

  // We can now make tests about the new tree
  v = {2, 2, 2};
  assert(f.covers_vector(idcs2, VType(std::move(v))));
  v = {7, 3, 1};
  assert(f.covers_vector(idcs2, VType(std::move(v))));
  v = {8, 3, 1};
  assert(not f.covers_vector(idcs2, VType(std::move(v))));
  v = {8, 5, 3};
  assert(not f.covers_vector(idcs2, VType(std::move(v))));

  // One more test so that covers needs to explore branches with shared
  // prefixes
  data = {{3, 2, 2}, {3, 4, 1}, {3, 2, 3}, {3, 4, 0}};
  auto idcs3 = f.add_vectors(std::move(vvtovv(data)));
  std::cout << f << std::endl;
  v = {3, 2, 3};
  assert(f.covers_vector(idcs3, VType(std::move(v))));
  v = {3, 2, 4};
  assert(not f.covers_vector(idcs3, VType(std::move(v))));
  v = {3, 4, 1};
  assert(f.covers_vector(idcs3, VType(std::move(v))));
  v = {3, 4, 2};
  assert(not f.covers_vector(idcs3, VType(std::move(v))));

  // Test that all children are ordered descending
  assert(f.check_child_order());

  // Test that the root of a (single-vector) ST is simulated by the
  // root of the original ST
  data = {{3, 4, 1}};
  auto idcs4 = f.add_vectors(std::move(vvtovv(data)));
  assert(f.check_simulation(idcs, idcs4));

  // Test that the root of a second (single-vector) ST is not simulated
  // by the root of the original ST
  data = {{4, 9, 9}};
  auto idcs5 = f.add_vectors(std::move(vvtovv(data)));
  assert(not f.check_simulation(idcs, idcs5));

  // Union tests: Check for vectors that are only in one ST
  // and a vector that is in none
  auto uRoot = f.st_union(idcs, idcs2);
  v = {7, 4, 1};
  assert(f.covers_vector(uRoot, VType(std::move(v))));
  v = {5, 5, 3};
  assert(f.covers_vector(uRoot, VType(std::move(v))));
  v = {8, 3, 1};
  assert(not f.covers_vector(uRoot, VType(std::move(v))));

  // Intersection test: Exact opposite of union
  auto iRoot = f.st_intersect(idcs, idcs2);
  v = {7, 4, 1};
  assert(not f.covers_vector(iRoot, VType(std::move(v))));
  v = {5, 5, 3};
  assert(not f.covers_vector(iRoot, VType(std::move(v))));
  v = {3, 5, 4};
  assert(f.covers_vector(iRoot, VType(std::move(v))));

  return 0;
}
