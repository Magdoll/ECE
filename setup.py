from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("findECE", ["findECE.c"])]
#ext_modules = [Extension("findECE", ["findECE.pyx"])]

setup(
		name = 'findECE',
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
)

