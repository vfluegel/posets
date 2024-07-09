downsets_dep = declare_dependency(include_directories: '.', dependencies : [boost_dep])

header_files = [
  'downsets/downsets/full_set.hh',
  'downsets/downsets/kdtree_backed.hh',
  'downsets/downsets/set_backed.hh',
  'downsets/downsets/vector_backed.hh',
  'downsets/downsets/vector_backed_bin.hh',
  'downsets/downsets/vector_backed_one_dim_split.hh',
  'downsets/downsets/vector_backed_one_dim_split_intersection_only.hh',
  'downsets/downsets/vector_or_kdtree_backed.hh',
  'downsets/downsets.hh',
  'downsets/utils/kdtree.hh',
  'downsets/utils/ref_ptr_cmp.hh',
  'downsets/utils/simd_traits.hh',
  'downsets/utils/vector_mm.hh',
  'downsets/vectors/array_backed.hh',
  'downsets/vectors/array_backed_sum.hh',
  'downsets/vectors/array_ptr_backed.hh',
  'downsets/vectors/simd_array_backed.hh',
  'downsets/vectors/simd_array_backed_sum.hh',
  'downsets/vectors/simd_array_ptr_backed.hh',
  'downsets/vectors/simd_po_res.hh',
  'downsets/vectors/simd_po_res_sum.hh',
  'downsets/vectors/simd_vector_backed.hh',
  'downsets/vectors/vector_backed.hh',
  'downsets/vectors.hh'
]

install_headers(header_files, preserve_path: true)
