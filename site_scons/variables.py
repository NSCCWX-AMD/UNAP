# -*- coding: utf-8 -*-
"""\
program Build-time Variables
----------------------------

This module defines the SCons variables used to control the compilation process
for program.

"""

import os
from SCons.Environment import Environment
import getpass
from SCons.Variables import (Variables, EnumVariable, PathVariable,
                             BoolVariable)
from SCons.Script import ARGUMENTS
import utils

### General Variables
SCONS_SITE_DIR = os.path.dirname(__file__)
PROJECT_DIR = os.path.dirname(SCONS_SITE_DIR)
PROJECT_BASE = os.path.dirname(PROJECT_DIR)
PROJECT_EXT_DIR = os.path.join(PROJECT_DIR, 'external')

program_vars = Variables('build_config.py', ARGUMENTS)
program_vars.AddVariables(
    # Project specific variables
    PathVariable('PROJECT_DIR', 'Path to the project directory', PROJECT_DIR,
                 PathVariable.PathIsDir),
    EnumVariable('PLATFORM',
                 'build platforms',
                 utils.ostype(),
                 allowed_values=('windows', 'linux', 'sw')),
    EnumVariable('BUILD_TYPE',
                 'Type of build',
                 'Opt',
                 allowed_values=('Opt', 'Debug', 'Prof')),
    EnumVariable('BUILD_ARCH',
                 'Build architecture',
                 '64',
                 allowed_values=('64', '32')),
    EnumVariable('PRECISION',
                 'Single/Double precision',
                 'DP',
                 allowed_values=('DP', 'SP')),
    EnumVariable('INT_TYPE', 'Integer size', '32',
                 allowed_values=('32', '64')),
    EnumVariable('FLOAT_TYPE', 'Float size', '64',
                 allowed_values=('32', '64')),
    BoolVariable('OMP', 'Use OpenMP multi-threading', False),
    EnumVariable('LIB_TYPE',
                 'library building type',
                 'static',
                 allowed_values=('shared', 'static', 'object')),
)

ostype = Environment(variables=program_vars)['PLATFORM']
print(ostype)

if ostype == "windows":
    print('ok windows')
    program_vars.AddVariables(
        ('CC', 'C compiler', 'gcc'),
        ('CXX', 'C++ compiler', 'g++'),
        ('F90', 'fortran90 compiler', 'gfortran'),
        ('CXX_LINKER', 'C++ linker', 'g++'),
        ('F_LINKER', 'fortran linker', 'gfortran'),
        ('MPI_LIB_NAME', 'MPI library name', 'mpi'),
        PathVariable('MPI_INC_PATH', 'Path to MPI headers',
                     'C:\Program Files\MPICH2\include',
                     PathVariable.PathIsDir),
        PathVariable('MPI_LIB_PATH', 'Path to MPI libraries',
                     'C:\Program Files\MPICH2\lib', PathVariable.PathIsDir),
    )
elif ostype == "linux":
    print('ok linux')
    program_vars.AddVariables(
        ('CC', 'C compiler', 'icc'),
        ('CXX', 'C++ compiler', 'icpc'),
        ('F90', 'fortran90 compiler', 'ifort'),
        ('CXX_LINKER', 'C++ linker', 'mpiicpc'),
        ('F_LINKER', 'fortran linker', 'mpiifort'),
        ('MPI_LIB_NAME', 'MPI library name', 'mpi'),
        PathVariable('MPI_INC_PATH', 'Path to MPI headers',
                     '/usr/sw-cluster/mpi2/include', PathVariable.PathIsDir),
        PathVariable('MPI_LIB_PATH', 'Path to MPI libraries',
                     '/usr/sw-cluster/mpi2/lib', PathVariable.PathIsDir),
    )
elif ostype == "sw":
    print('ok sw')
    program_vars.AddVariables(
        ('CC_HOST', 'c compiler on host', 'sw5cc'),
        ('CC_SLAVE', 'c compiler on slave', 'sw5cc'),
        ('CC', 'C compiler', 'swgcc'),
        ('CXX', 'C++ compiler', 'swg++453'),
        ('F90', 'fortran90 compiler', 'mpif90'),
        ('CXX_LINKER', 'C++ linker', 'swld453'),
        ('F_LINKER', 'fortran linker', 'swld453-fort'),
        ('MPI_LIB_NAME', 'MPI library name', 'mpi'),
        PathVariable('MPI_INC_PATH', 'Path to MPI headers',
                     '/usr/sw-mpp/mpi2/include', PathVariable.PathIsDir),
        PathVariable('MPI_LIB_PATH', 'Path to MPI libraries',
                     '/usr/sw-mpp/mpi2/lib', PathVariable.PathIsDir),
        BoolVariable('ATHREAD', 'Use Shenwei multi-threading', False),
    )
else:
    print('Unknown ostype')


def init_dependent_vars(env):
    """Initialize dependent variables based on user configuration"""
    from SCons.Script import Mkdir
    # ostype = utils.ostype()
    ostype = env['PLATFORM']
    prj_dir = env['PROJECT_DIR']

    BUILD_OPTION = (ostype + env['BUILD_ARCH'] + env['CXX'] +
                    env['PRECISION'] + env['BUILD_TYPE'])
    PLATFORM_INSTALL = os.path.join(prj_dir, 'install', BUILD_OPTION)
    BIN_PLATFORM_INSTALL = os.path.join(PLATFORM_INSTALL, 'bin')
    LIB_PLATFORM_INSTALL = os.path.join(PLATFORM_INSTALL, 'lib')
    INC_PLATFORM_INSTALL = os.path.join(PLATFORM_INSTALL, 'include')
    LIB_SRC = os.path.join(prj_dir, 'src')
    # PROJECT_INC_DIR = os.path.join(prj_dir, 'build', 'lnInclude')
    PROJECT_INC_DIR = INC_PLATFORM_INSTALL
    EXTERNAL_DIR = os.path.join(prj_dir, 'external')
    EXTERNAL_WINDOWS_DIR = os.path.join(EXTERNAL_DIR, 'windows')
    env.Append(
        BUILD_OPTION=BUILD_OPTION,
        BIN_PLATFORM_INSTALL=BIN_PLATFORM_INSTALL,
        LIB_PLATFORM_INSTALL=LIB_PLATFORM_INSTALL,
        INC_PLATFORM_INSTALL=INC_PLATFORM_INSTALL,
        LIB_SRC=LIB_SRC,
        PROJECT_INC_DIR=PROJECT_INC_DIR,
        EXTERNAL_DIR=EXTERNAL_DIR,
        EXTERNAL_WINDOWS_DIR=EXTERNAL_WINDOWS_DIR,
    )
