# -*- coding: utf-8 -*-
"""\
SCons build rules
"""

import os

cxx_source_files = []
c_source_files = []
fortran_source_files = []
chost_source_files = []
cslave_source_files = []
cxxhost_source_files = []


def build_object(baseenv,
                 sources,
                 program_inc,
                 program_libs,
                 sources_type='none',
                 prepend_args=None,
                 append_args=None):
    libenv = baseenv.Clone()
    libenv.Prepend(CPPPATH=program_inc)

    if sources_type == 'host':
        libenv.Replace(
            CCCOM='$CC_HOST -host -OPT:IEEE_arith=2 -c -o $TARGET $SOURCES')
    elif sources_type == 'slave':
        libenv.Replace(
            CCCOM=
            '$CC_SLAVE -slave -OPT:IEEE_arith=2 -msimd -c -o $TARGET $SOURCES')

    objs = libenv.Object(source=sources)
    return objs


def build_lib(baseenv,
              target,
              sources,
              program_inc,
              program_libs,
              prepend_args=None,
              append_args=None):
    """Build a shared library

    Args:
        baseenv (env): program SCons build environment
        target (str): Name of the build target
        sources (list): List of sources for this target
        program_inc (list): List of program include paths
        program_libs (list): List of libraries to be linked

        prepend_args (dict): Set of (key, value) pairs to be prepended
        append_args (dict): Set of (key, value) pairs to be appended
    """
    libenv = baseenv.Clone()
    lib_type = libenv['LIB_TYPE']
    lib_src = libenv['LIB_SRC']
    inc_dirs = [os.path.join(lib_src, d) for d in program_inc]
    libenv.Prepend(CPPPATH=inc_dirs)
    libenv.Append(LIBS=program_libs)
    libenv.Append(LIBPATH=libenv['LIBPATH_COMMON'] + libenv['LIBPATH_LIBS'])

    if lib_type == "shared":
        exe = libenv.SharedLibrary(target=target, source=sources)
    elif lib_type == "static":
        exe = libenv.StaticLibrary(target=target, source=sources)

    install_dir = libenv['LIB_PLATFORM_INSTALL']
    libenv.Alias('install', install_dir)
    libenv.Install(install_dir, exe)

    if prepend_args is not None:
        libenv.Prepend(**prepend_args)
    if append_args is not None:
        libenv.Append(**append_args)

    return libenv


def build_app(baseenv,
              target,
              sources,
              program_inc,
              program_libs,
              prepend_args=None,
              append_args=None,
              linker=None):
    """Build an executable application

    Args:
        baseenv (env): program SCons build environment
        target (str): Name of the build target
        sources (list): List of sources for this target
        program_inc (list): List of program include paths
        program_libs (list): List of libraries to be linked

        prepend_args (dict): Set of (key, value) pairs to be prepended
        append_args (dict): Set of (key, value) pairs to be appended
    """
    appenv = baseenv.Clone()

    lib_src = appenv['LIB_SRC']
    inc_dirs = [os.path.join(lib_src, d) for d in program_inc]
    appenv.Prepend(CPPPATH=inc_dirs)
    appenv.Prepend(F90PATH=inc_dirs)
    appenv.Append(LIBS=program_libs)
    appenv.Append(LIBPATH=appenv['LIBPATH_COMMON'] + appenv['LIBPATH_APPS'])

    exe = appenv.Program(target=target, source=sources)
    install_dir = appenv['BIN_PLATFORM_INSTALL']
    appenv.Alias('install', install_dir)
    appenv.Install(install_dir, exe)

    if prepend_args is not None:
        appenv.Prepend(**prepend_args)
    if append_args is not None:
        appenv.Append(**append_args)

    if linker is not None:
        appenv['LINK'] = linker

    return appenv


def build_lninclude(env):
    """Create lnInclude directories

    Args:
        env (Environment): program SCons build environment
    """

    from SCons.Script import Copy, Mkdir, Dir

    ostype = "windows" if env['PLATFORM'] == "windows" else "posix"
    include_patterns = [".hpp", ".H", ".hxx", ".h", ".hh"]

    inc_env = env.Clone()
    inc_dir = inc_env['PROJECT_INC_DIR']
    Mkdir(inc_dir)
    src_dir = inc_env['LIB_SRC']
    if os.path.exists(inc_dir):
        inc_env.Alias("install", inc_dir)

    for root, dlist, files in os.walk(src_dir):
        if "lnInclude" in dlist:
            dlist.remove("lnInclude")

        dbase = os.path.basename(root)
        if "OSspecific" in dbase:
            for d in dlist:
                if d != ostype:
                    dlist.remove(d)

        for f in files:
            if any(f.endswith(pat) for pat in include_patterns):
                src = os.path.join(root, f)
                dest = os.path.join(inc_dir, f)
                yield inc_env.Install(inc_dir, src)


def add_source_files(source_files, all_source_files):
    for filename in source_files:
        fullDir = os.path.join(os.getcwd(), filename)
        all_source_files.append(fullDir)
