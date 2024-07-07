project ('Downset Library', ['cpp'],
         default_options : ['cpp_std=c++20', 'b_ndebug=if-release'],
         license : 'GPL-3.0')

cpp = meson.get_compiler('cpp')

if cpp.get_id() == 'clang'
  extra_args = []
else
  extra_args = ['-ftrack-macro-expansion=0',
                '-fno-diagnostics-show-caret']
endif

add_global_arguments (extra_args, language: 'cpp')

boost_dep = dependency ('boost')

if cpp.has_header ('experimental/simd')
  stdsimd_dep = declare_dependency ()
else
  stdsimd_dep = dependency ('stdsimd')
endif

subdir ('include')
subdir ('tests')
