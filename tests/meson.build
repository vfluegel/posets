tests_exe = executable ('tests', 'tests.cc',
                        dependencies : [posets_dep ],
                        c_args : '-O0' )  # Deactivate optimization, otherwise it takes ages.

kdtests_exe = executable ('kdtests', 'kdtests.cc',
                          dependencies : [posets_dep ])

sttests_exe = executable ('sttests', 'sttests.cc',
                          dependencies : [posets_dep ])

simdbm_exe = executable ('simd-bm', 'simd-bm.cc',
                         dependencies : [posets_dep ])

downsetbm_exe = executable ('downset-bm', 'downset-bm.cc',
                            dependencies : [posets_dep ])

test('antichains/vectors random', downsetbm_exe,
     args : ['all', 'all', '--params=build=20,query=30,transfer=100,intersection=10,union=20' ])
test('antichains/vectors implementations', tests_exe, args : ['all', 'all'])
test('antichains/vectors kdtests', kdtests_exe)
test('antichains/vectors sttests', sttests_exe)
