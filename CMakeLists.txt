# ======================================================================
#  larana main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larana/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(larana VERSION 09.04.02 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED PUBLIC)
find_package(art_root_io REQUIRED PUBLIC)
find_package(nug4 REQUIRED PUBLIC)
find_package(nurandom REQUIRED PUBLIC)
find_package(nutools REQUIRED PUBLIC)
find_package(PostgreSQL REQUIRED PUBLIC)
find_package (Eigen3 3.3...3.3.9 REQUIRED NO_MODULE)
find_package(ROOT COMPONENTS Core RIO Physics Tree Hist MathCore TMVA REQUIRED PUBLIC)
find_package(CLHEP COMPONENTS Random REQUIRED PUBLIC)

find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( larreco )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

include_directories( $ENV{EIGEN_INC} )

# source
add_subdirectory(larana)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

