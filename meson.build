project ('Posets Library', ['cpp'],
         default_options : ['cpp_std=c++20', 'b_ndebug=if-release'],
         license : 'GPL-3.0')

if not meson.is_subproject()
  cpp = meson.get_compiler('cpp')
  if cpp.get_id() == 'clang'
    extra_args = ['-fexperimental-library']
  else
    extra_args = ['-ftrack-macro-expansion=0',
                  '-fno-diagnostics-show-caret',
                  '-fconcepts-diagnostics-depth=2']
  endif
  add_global_arguments (extra_args, language: 'cpp')
endif

boost_dep = dependency ('boost')

subdir ('include')
subdir ('tests')
