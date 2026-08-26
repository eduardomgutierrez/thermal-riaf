project( 'adaf',
         ['cpp', 'c'],
         default_options : [ 'cpp_std=c++17',
                             'buildtype=release' ],
         version : '0.1' )

gsl_dep = dependency('gsl')
boost_dep = dependency('boost')
omp_dep = dependency('openmp')
# Request the C++ interface through Meson's HDF5 resolver. On Debian/Ubuntu,
# libhdf5_cpp lives in an HDF5-specific library directory and a bare
# compiler.find_library('hdf5_cpp') cannot find it. This dependency also links
# the underlying HDF5 C library.
hdf5_cpp_dep = dependency('hdf5', language : 'cpp')

subdir('src/lib/nrMath')
subdir('src/lib/fmath')
subdir('src/lib/fparameters')
subdir('src/lib/inout')
subdir('src/lib/fparticle')
subdir('src/lib/finjection')
subdir('src/lib/flosses')
subdir('src/lib/fluminosities')
subdir('src/adaf/')
