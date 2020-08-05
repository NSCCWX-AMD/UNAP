# -*- mode: python -*-

import os
import utils
from variables import program_vars, init_dependent_vars
from compiler import update_compiler_settings
from scons_color import color_prints

### Initialize toolsets based on operating system
ostype = Environment(variables = program_vars)['PLATFORM']
tools = ['default']
if ostype == 'windows':
    tools += ['mingw']

### Base SCons environment
env = Environment(variables = program_vars,
                  tools = tools,
                  ENV = os.environ)
Help(program_vars.GenerateHelpText(env))
init_dependent_vars(env)
update_compiler_settings(env)

### simpler compiling message
### if you want detailed information, comment next line
color_prints(env)

### Isolate build environments based on build options
build_dir = os.path.join(
    Dir("#").abspath,"build",env['BUILD_OPTION'])

program_src = [
    'src',
    'test'
]

for d in program_src:
    SConscript('%s/SConscript'%d,
               exports=['env'],
               src_dir=Dir("#").srcnode().abspath,
               variant_dir=build_dir)

### Remove buid directory when cleaning
Clean(".", build_dir)
