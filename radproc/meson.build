project( 'adaf',
         ['cpp', 'c'],
         default_options : [ 'cpp_std=c++17',
                             'buildtype=release' ],
         version : '0.1' )

gsl_dep = dependency('gsl')
boost_dep = dependency('boost')
omp_dep = dependency('openmp')

subdir('src/lib/nrMath')
subdir('src/lib/fmath')
subdir('src/lib/fparameters')
subdir('src/lib/inout')
subdir('src/lib/fparticle')
subdir('src/lib/finjection')
subdir('src/lib/flosses')
subdir('src/lib/fluminosities')
subdir('src/adaf/')
