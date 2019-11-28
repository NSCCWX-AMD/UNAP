# GLOVAR.PY
import os

# FUNCTION DEFINATION
# initialize variable
def _init():
	global global_files
	global_files = []

# add cppfiles
def add_files(cxxfile):
	for filename in cxxfile:
		fullDir = os.getcwd() + '/' + filename
		# relPath = os.path.relpath(fullDir, root_path)
		global_files.append(fullDir)

# return variable
def get_files(defValue=None):
	try:
		return global_files
	except KeyError:
		return defValue

# build symbolic link
def Symlink (target, source, env):
	os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))

# add target and source for symbolic link function
def get_source_headerfiles(current_path):

	source_headerfiles = []
	source_path       = []

	for (root, dirnames, filenames) in os.walk(current_path):
		for filename in filenames:
			if filename.endswith((".hpp",".h")):
				source_path = root + '/'
				source_path = source_path + filename
				source_headerfiles.append(source_path)
	return source_headerfiles

def get_target_headerfiles(source_header_files):

	target_headerfiles = []
	target_path        = []

	for i in range(len(source_header_files)):
		filename = os.path.basename(source_header_files[i])
		target_path = include_path + '/' + filename
		target_headerfiles.append(target_path)
	return target_headerfiles

#build symbolic link
def build_Include(env_local):
	source_header_files = get_source_headerfiles(source_path)
	target_header_files = get_target_headerfiles(source_header_files)

	for i in range(len(source_header_files)):
		env_local.Command(target_header_files[i], source_header_files[i], Symlink)
	return True

# PATH DEFINATION
fullName = ''
root_path          = os.getcwd()
include_path       = ''
source_path        = ''
build_path         = ''
platform           = ''

