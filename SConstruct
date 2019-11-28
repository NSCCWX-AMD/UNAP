import os
import GLOVAR as gl

# initialization
gl._init()

software_home = '/home/export/online1/amd_dev1/software'

# Set our required environment
libraries 		= ['mpi']
library_paths 	= [software_home + '/MPICH/lib']
cppDefines 		= {}
cppFlags 		= ['-fPIC', '-DSCALAR_FLOAT64', '-g', '-DSWTIMER']
debugFlags      = ['-g', '-O0', '-DDEBUG', '-Wall', '-Wextra', '-Werror']
optFlags		= ['-O3']
f90Flags		= ['-fPIC', '-g', '-cpp']

# define the attributes of the build environment
env = Environment(ENV   = os.environ)
env.Append(LIBS 	    = libraries)
env.Append(LIBPATH 		= library_paths)
env.Append(CPPDEFINES 	= cppDefines)
env.Append(CPPFLAGS 	= cppFlags)
env.Append(F90FLAGS 	= f90Flags)
env.Append(CCCOMSTR     = "CC $SOURCES")
#env.Append(LINKCOMSTR   = "LINK $TARGET")

# get debug flag, the default is -O3
debug = ARGUMENTS.get('debug', '')
DBGFLAG = ''

if debug == 'true':
	env.Append(CPPFLAGS = debugFlags)
	DBGFLAG = 'Debug'
else:
	env.Append(CPPFLAGS = optFlags)
	env.Append(F90FLAGS = optFlags)
	DBGFLAG = 'Opt'

SWOPT = ''
INTFLAG = '-DLABEL_INT64'
INTSize = ''

if INTFLAG == '-DLABEL_INT32':
    INTSize = 'Int32'
elif INTFLAG == '-DLABEL_INT64':
    INTSize = 'Int64'
else:
    print('Error: label size is not defined!')

env.Append(CPPFLAGS = INTFLAG)

# get platform information, the default is x86
# if sw preferred, using "scons pla=sw"
gl.platform = ARGUMENTS.get('pla', 'x86')

if gl.platform == 'sw':
	env['CC']   = 'sw5cc'
	env['CXX']  = 'swg++453'
	env['F90']  = 'mpif90'
	env['AR']   = 'swar'
	env['LINK'] = 'swg++453'
	env.Append(CPPPATH = ['/usr/sw-mpp/mpi2/include'])
else:
	env['CC']   = 'icc'
	env['CXX']  = 'icpc'
	env['F90']  = 'ifort'
	env['AR']   = 'ar'
	env.Append(CPPPATH = [software_home + '/MPICH/include'])
	env.Append(F90PATH = [software_home + '/MPICH/include'])

# env['LINKFLAGS'] = '-lstdc++'

Export('env')

src_name = 'src'
build_name = 'build'
test_name = 'test'

gl.include_path = gl.root_path + '/' + build_name + '/Include'
print(gl.root_path)

# add include path to -I
env.Append(CPPPATH = gl.include_path)

CXXNAME = os.path.basename(env['CXX'])

gl.fullName = CXXNAME + DBGFLAG + INTSize + SWOPT
gl.build_path = gl.root_path + '/' + build_name + '/' + gl.fullName
gl.source_path = gl.root_path + '/' + src_name
build_src_path = gl.build_path + '/' + src_name
build_test_path = gl.build_path + '/' + test_name

# out of source compiling settings
env.VariantDir(build_src_path, gl.source_path, duplicate=1)
env.SConscript(build_src_path + '/SConscript', duplicate=0, exports = 'env' )

env.VariantDir(build_test_path, gl.root_path + '/' + test_name, duplicate=1)
env.SConscript(build_test_path + '/SConscript', duplicate=0, exports = 'env' )
