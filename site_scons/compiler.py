# -*- coding: utf-8 -*-
"""\
Compiler settings
-------------------------
"""

import os

generalflags = dict(
    general="-fPIC",
    # warnings = "-Wall -Wextra",
    warnings="",
    debug="-O0 -ggdb3 -DDEBUG",
    prof="-O2 -pg",
    opt="-O3")

gcc_flags = dict(**generalflags)

intel_flags = dict(**generalflags)
intel_flags['warnings'] = "-wd327,654,819,1125,1476,1505,1572"

clang_flags = dict(**generalflags)
# clang_flags['warnings'] = "-Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wno-overloaded-virtual -Wno-unused-comparison -Wno-deprecated-register"

#gcc_flags['warnings'] = clang_flags['warnings']

_compiler_map = {
    'gfortran': 'g++',
    'gcc': 'g++',
    'g++': 'g++',
    'cc': 'c++',
    'c++': 'c++',
    'clang': 'clang++',
    'clang++': 'clang++',
    'icc': 'icpc',
    'intel': 'icpc',
    'icpc': 'icpc',
    'swgcc': 'g++'
}

_compiler_flags_map = {
    'g++': gcc_flags,
    'icpc': intel_flags,
    'clang++': clang_flags,
    'c++': clang_flags
}


def windows_flags(env):
    """Windows specific compiler flags"""
    env.Append(
        CPPPATH=[
            env['MPI_INC_PATH'],
        ],
        F90PATH=[env['MPI_INC_PATH']],
        LIBPATH=[
            env['MPI_LIB_PATH'],
        ],
        # LIBPATH_COMMON=[env['MPI_LIB_PATH']]
    )
    # env.Prepend(LINKFLAGS='-Xlinker --add-needed')
    # env.Prepend(LINKFLAGS='-Xlinker --no-as-needed')
    # env.Prepend(LINKFLAGS='-Wl,--whole-archive')
    # env.Prepend(LINKFLAGS='-Wl,-no-whole-archive')


def linux_flags(env):
    """Linux specific compiler flags"""
    env.Append(LIBPATH_COMMON=[env['MPI_LIB_PATH']],
               CPPPATH=[env['MPI_INC_PATH']],
               F90PATH=[env['MPI_INC_PATH']])

    # env.Prepend(LINKFLAGS = '-Xlinker --add-needed')
    # env.Prepend(LINKFLAGS = '-Xlinker --no-as-needed')


def darwin_flags(env):
    """Darwin specific compiler flags"""
    env.Append(CXXFLAGS='-std=c++11')

    env.Append(LIBPATH_COMMON=[env['MPI_LIB_PATH']],
               CPPPATH=[env['MPI_INC_PATH']],
               LIBS=[
                   'pthread',
               ])

    env.Prepend(LINKFLAGS='-undefined dynamic_lookup')


def sunway_flags(env):
    """Sunway specific compiler flags"""
    print('sunway has nothing special')
    env.Append(LIBPATH_COMMON=[env['MPI_LIB_PATH']],
               CPPPATH=[env['MPI_INC_PATH']],
               F90PATH=[env['MPI_INC_PATH']])

    # env.Prepend(LINKFLAGS = '-Xlinker --add-needed')
    # env.Prepend(LINKFLAGS = '-Xlinker --no-as-needed')


_arch_map = dict(
    windows=windows_flags,
    linux=linux_flags,
    darwin=darwin_flags,
    sw=sunway_flags,
)


def update_compiler_settings(env):
    """Update compiler flags and other settings

    Args:
        env (SCons.Environment): Environment to modify
    """
    ostype = env['PLATFORM']
    c_compiler = env['CC']
    cxx_compiler = env['CXX']
    fort_compiler = env['F90']
    btype = env['BUILD_TYPE']
    flist = ['general', 'warnings', btype.lower()]
    keys = _compiler_map.keys()
    cname = os.path.basename(c_compiler)
    for k in keys:
        if k in cname:
            cname = k
            break
    cflags = _compiler_flags_map[_compiler_map[cname]]
    # Compiler flags
    env.Append(CCFLAGS='-DLABEL_INT' + env['INT_TYPE'])
    env.Append(CCFLAGS='-DSCALAR_FLOAT' + env['FLOAT_TYPE'])

    env.Append(CCFLAGS='-D' + ostype)
    env.Append(CCFLAGS='-DWM_' + env['PRECISION'])
    for k in flist:
        env.Append(CCFLAGS=cflags[k].split())

    if env['PLATFORM'] == 'sw':
        if env['ATHREAD']:
            env.Append(CCFLAGS='-DSW_SLAVE')

    # Linker flags
    env.Append(LIBPATH_COMMON=[env['LIB_PLATFORM_INSTALL']],
               LIBPATH_APPS=[],
               LIBPATH_LIBS=[])
    # env.Append(RPATH = [env['LIB_PLATFORM_INSTALL']])

    # Global include directory
    env.Append(CPPPATH=[env['PROJECT_INC_DIR']])
    env.Append(F90PATH=[env['PROJECT_INC_DIR']])

    arch_func = _arch_map.get(ostype, None)
    if arch_func is not None:
        arch_func(env)

    # LEX replacement
    # env.Replace(LEXCOM='$LEX $LEXFLAGS -o$TARGET -f $SOURCES')

    # Fix MPI flags
    mpi_lib = env['MPI_LIB_NAME']
    env.Append(LIBS=mpi_lib)
    if "mpich" in mpi_lib.lower():
        env.Append(CCFLAGS='-DMPICH_SKIP_MPICXX')
    else:
        env.Append(CCFLAGS='-DOMPI_SKIP_MPICXX')

    env['F90FLAGS'] = env['CCFLAGS']

    # make gfortran support preprocessor
    env.Append(F90FLAGS='-cpp')
