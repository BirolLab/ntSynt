project('ntSynt', 'cpp',
        version : '1.0.5',
        license : 'GPL-3',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true' ])

compiler = meson.get_compiler('cpp')
btllib_dep = compiler.find_library('btllib')

openmp_dep = dependency('openmp')

include_dirs = [
    include_directories('src/vendor/argparse/include'),
]

rootpath = meson.source_root()

src_files = [ 'ntsynt_make_common_bf.cpp' ]
deps = [ btllib_dep, openmp_dep ]

link_args = [ '-lbtllib' ]

clang_format = find_program('clang-format', required : false)
if clang_format.found()
  run_target('format', command: join_paths(rootpath, 'src/scripts/format'))
endif

# Source files and scripts
subdir('src')
subdir('bin')

# Install required files from ntJoin submodule
meson.add_install_script('scripts/install-ntjoin')