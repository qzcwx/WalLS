from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_module1 = Extension(
    "demo",
    ["demo.pyx"],
    extra_compile_args=['-pipe'],
#    extra_compile_args=['-fopenmp'],
#    extra_link_args=['-fopenmp'],
#    language="c++",
)

ext_module2 = Extension(
    "LocalSearch",
    ["LocalSearch.pyx"],
#    extra_compile_args=['-fopenmp'],
#    extra_link_args=['-fopenmp'],
    extra_compile_args=['-pipe'],    
    language="c++",
)


ext_module3 = Extension(
    "nkLandscape",
    ["nkLandscape.pyx"],
#    extra_compile_args=['-fopenmp'],
#    extra_link_args=['-fopenmp'],
   language="c++",
)

ext_module4 = Extension(
    "nkqLandscape",
    ["nkqLandscape.pyx"],
#    extra_compile_args=['-fopenmp'],
#    extra_link_args=['-fopenmp'],
   language="c++",
)

ext_module5 = Extension(
    "tool",
    ["tool.pyx"],
#    extra_compile_args=['-fopenmp'],
#    extra_link_args=['-fopenmp'],
)

ext_module6 = Extension(
    "WalshAnalysis",
    ["WalshAnalysis.pyx"],
#    extra_compile_args=['-fopenmp'],
#    extra_link_args=['-fopenmp'],
#    language="c++",
)

ext_module7 = Extension(
    "individual",
    ["individual.pyx"],
#    extra_compile_args=['-fopenmp'],
#    extra_link_args=['-fopenmp'],
    language="c++",
)

ext_module8 = Extension(
    "LocalOptima",
    ["LocalOptima.pyx"],
    language="c++",
)

setup(
    name = 'Hello world app',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_module1, ext_module2, ext_module3, ext_module4, ext_module5, ext_module6, ext_module7, ext_module8],
)
