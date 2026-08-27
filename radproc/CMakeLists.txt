cmake_minimum_required(VERSION 3.16)
project(adaf VERSION 0.1 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
set(LIB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/src/lib)

add_library(nrMath STATIC
    src/lib/nrMath/bessel.cpp src/lib/nrMath/bisection.cpp
    src/lib/nrMath/brent.cpp src/lib/nrMath/interpol.cpp
    src/lib/nrMath/laguerre.cpp src/lib/nrMath/broydn.cpp
    src/lib/nrMath/erff.cpp src/lib/nrMath/integrators.cpp
    src/lib/nrMath/lnsrch.cpp src/lib/nrMath/lubksb.cpp
    src/lib/nrMath/ludcmp.cpp src/lib/nrMath/newt.cpp
    src/lib/nrMath/nrutil.cpp src/lib/nrMath/ran1.cpp
    src/lib/nrMath/zbrac.cpp src/lib/nrMath/zbrent.cpp)
target_include_directories(nrMath PUBLIC ${LIB_ROOT}/nrMath ${LIB_ROOT})
target_link_libraries(nrMath PUBLIC GSL::gsl GSL::gslcblas)

add_library(fmath STATIC
    src/lib/fmath/fbisection.cpp src/lib/fmath/matrixInit.cpp
    src/lib/fmath/tridiagonal.cpp src/lib/fmath/interpolation.cpp
    src/lib/fmath/mathFunctions.cpp src/lib/fmath/RungeKutta.cpp)
target_include_directories(fmath PUBLIC ${LIB_ROOT}/fmath ${LIB_ROOT})

add_library(fparameters STATIC
    src/lib/fparameters/Dimension.cpp src/lib/fparameters/DimensionIterator.cpp
    src/lib/fparameters/parameters.cpp src/lib/fparameters/ParamSpace.cpp
    src/lib/fparameters/ParamSpaceValues.cpp src/lib/fparameters/SpaceCoord.cpp
    src/lib/fparameters/SpaceIterator.cpp)
target_include_directories(fparameters PUBLIC ${LIB_ROOT}/fparameters ${LIB_ROOT})
target_link_libraries(fparameters PUBLIC fmath Boost::boost)

add_library(inout STATIC src/lib/inout/ioutil.cpp)
target_include_directories(inout PUBLIC ${LIB_ROOT}/inout ${LIB_ROOT})
target_link_libraries(inout PUBLIC fparameters Boost::boost)

add_library(fparticle STATIC src/lib/fparticle/Particle.cpp)
target_include_directories(fparticle PUBLIC ${LIB_ROOT}/fparticle ${LIB_ROOT})
target_link_libraries(fparticle PUBLIC fparameters Boost::boost)

add_library(finjection STATIC
    src/lib/finjection/muonInj.cpp src/lib/finjection/neutrinoInj.cpp
    src/lib/finjection/neutronDecay.cpp src/lib/finjection/neutronPp.cpp
    src/lib/finjection/pairBH.cpp src/lib/finjection/pairMuonDecay.cpp
    src/lib/finjection/ppPionInj.cpp src/lib/finjection/pairgg.cpp
    src/lib/finjection/pgammaPionInj.cpp src/lib/finjection/neutronPgamma.cpp
    src/lib/finjection/pairInjectionExact.cpp)
target_include_directories(finjection PUBLIC ${LIB_ROOT}/finjection ${LIB_ROOT})
target_link_libraries(finjection PUBLIC Boost::boost)

add_library(flosses STATIC
    src/lib/flosses/crossSectionInel.cpp src/lib/flosses/lossesBrem.cpp
    src/lib/flosses/lossesHadronics.cpp src/lib/flosses/lossesIC.cpp
    src/lib/flosses/lossesPhotoHadronic.cpp src/lib/flosses/lossesSyn.cpp
    src/lib/flosses/nonThermalLosses.cpp)
target_include_directories(flosses PUBLIC ${LIB_ROOT}/flosses ${LIB_ROOT})

add_library(fluminosities STATIC
    src/lib/fluminosities/blackBody.cpp src/lib/fluminosities/compton.cpp
    src/lib/fluminosities/luminosityHadronic.cpp
    src/lib/fluminosities/luminosityIC.cpp
    src/lib/fluminosities/luminosityNTHadronic.cpp
    src/lib/fluminosities/luminosityPairAnnihilation.cpp
    src/lib/fluminosities/luminosityPhotoHadronic.cpp
    src/lib/fluminosities/luminositySynchrotron.cpp
    src/lib/fluminosities/opticalDepthSSA.cpp
    src/lib/fluminosities/probaprox.cpp src/lib/fluminosities/probexact.cpp
    src/lib/fluminosities/reflection.cpp
    src/lib/fluminosities/thermalBremss.cpp src/lib/fluminosities/thermalSync.cpp)
target_include_directories(fluminosities PUBLIC ${LIB_ROOT}/fluminosities ${LIB_ROOT})
target_link_libraries(fluminosities PUBLIC fparticle nrMath fmath flosses)

add_executable(adaf
    src/adaf/modelParameters.cpp src/adaf/thermalProcesses.cpp
    src/adaf/secondariesProcesses.cpp src/adaf/absorption.cpp
    src/adaf/adafFunctions.cpp src/adaf/redshiftFunction.cpp
    src/adaf/NTdistribution.cpp src/adaf/thermalCompton.cpp
    src/adaf/read.cpp src/adaf/State.cpp src/adaf/adafParameters.cpp
    src/adaf/NTtimescales.cpp src/adaf/NTinjection.cpp
    src/adaf/messages.cpp src/adaf/write.cpp src/adaf/thermalDistribution.cpp
    src/adaf/main.cpp src/adaf/losses.cpp src/adaf/comptonScattMatrix.cpp
    src/adaf/NTradiation.cpp src/adaf/maxwellJuttner.cpp)
target_include_directories(adaf PRIVATE src/adaf ${LIB_ROOT})
target_link_libraries(adaf PRIVATE
    nrMath fmath fparameters inout fparticle finjection flosses fluminosities
    OpenMP::OpenMP_CXX GSL::gsl GSL::gslcblas Boost::boost)
target_compile_options(adaf PRIVATE -O2 -march=native -ffast-math)
set_target_properties(adaf PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src/adaf")
