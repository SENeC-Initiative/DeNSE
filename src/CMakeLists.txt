cmake_minimum_required( VERSION 2.6 )
set( CMAKE_LEGACY_CYGWIN_WIN32 0 ) # for cygwin

# add cmake modules: for all `include(...)` first look here
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

if ( COMMAND cmake_policy )

#  cmake_policy( SET CMP0003 NEW )
cmake_policy( SET CMP0012 NEW ) # accept "TRUE" "FALSE" in if statements

endif()

set( SUCCESS 1 )


###############################################################################
#################                   Project                   #################
###############################################################################

project (cgrowth)

set( PY_NAME NetGrowth )

# Version number.
set (cgrowth_VERSION_MAJOR 0)
set (cgrowth_VERSION_MINOR 1)
set( cgrowth_VERSION_PATCHLEVEL 0 )

set( cgrowth_VERSION_VERSION "${cgrowth_VERSION_MAJOR}.${cgrowth_VERSION_MINOR}.${cgrowth_VERSION_PATCHLEVEL}" )


###############################################################################
#################            User-defined options             #################
###############################################################################

set( with-geos ON CACHE STRING "Give directory to CGAL library. [default=ON]" )
set( with-mpi OFF CACHE STRING "Request compilation with MPI. Optionally give directory with MPI installation. [default=OFF]" )
set( with-openmp ON CACHE STRING "Enable OpenMP multithreading. Optional: set OMP flag. [default=ON]" )
set( with-python ON CACHE STRING "Build python bindings. To set a specific Python, set install path. [default=ON]" )
set( CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Don't override flags")
option( cythonize-pybindings "Use Cython to cythonize `_growth.pyx`. If OFF, PyNEST has to be build from a pre-cythonized `_growth.pyx`. [default=ON]" ON )
option( with-docs OFF "Compile documentation (requires `doxygen` and `sphinx`. [default=OFF]" )
option( with-debug OFF "Compile with debug symbols. [default=OFF]" )
option( python-develop OFF "Python compilere allows to modify the code. [default=OFF]" )

# These includes publish function names.
include( ConfigureSummary )
include( GetTriple )
include( ProcessOptions )

# get triples arch-vendor-os
get_host_triple( CGROWTH_HOST_TRIPLE CGROWTH_HOST_ARCH CGROWTH_HOST_VENDOR CGROWTH_HOST_OS )
get_target_triple( CGROWTH_TARGET_TRIPLE CGROWTH_TARGET_ARCH CGROWTH_TARGET_VENDOR CGROWTH_TARGET_OS )

# process options
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
cgrowth_process_with_optimize()
cgrowth_process_with_debug()
cgrowth_process_with_warning()
cgrowth_process_with_libraries()
cgrowth_process_with_includes()
cgrowth_process_with_defines()
cgrowth_process_with_python()
cgrowth_process_with_geos()
cgrowth_process_with_openmp()
cgrowth_process_with_mpi()
cgrowth_process_with_docs()
cgrowth_get_color_flags()
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CGROWTH_C_COLOR_FLAGS}" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGROWTH_CXX_COLOR_FLAGS}" )

# check installation directory
message( "${CMAKE_INSTALL_PREFIX}" )
if ( "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" )
  # no prefix given?
  message( FATAL_ERROR "Do not install into '${CMAKE_INSTALL_PREFIX}'. "
                    "Please set -DCMAKE_INSTALL_PREFIX:PATH=<install prefix>." )
endif ()
if (NOT DEFINED CMAKE_INSTALL_LIBDIR)
  if ("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
    set(CMAKE_INSTALL_LIBDIR "lib64")
  else ()
    set(CMAKE_INSTALL_LIBDIR "lib")
  endif ()
endif ()
if (NOT DEFINED CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "bin")
endif ()


###############################################################################
#################        Define Subdirectories here           #################
###############################################################################

set (SOURCE_FILES "module.hpp" "module.cpp")
add_subdirectory( elements )
add_subdirectory( kernel )
add_subdirectory( libgrowth )
add_subdirectory( io )
add_subdirectory( models )
add_subdirectory( spatial )
if ( HAVE_PYTHON )
  add_subdirectory( pymodule )
endif ()
if ( WITH_DOCS )
  add_subdirectory( docs )
endif ()


###############################################################################
#################             Summary of flags                #################
###############################################################################

# all compiler flags
if ( NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "None" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS}" )
elseif ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}   ${CMAKE_C_FLAGS_DEBUG}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif ( ${CMAKE_BUILD_TYPE} STREQUAL "Release" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}   ${CMAKE_C_FLAGS_RELEASE}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}" )
elseif ( ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}   ${CMAKE_C_FLAGS_RELWITHDEBINFO}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )
elseif ( ${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}   ${CMAKE_C_FLAGS_MINSIZEREL}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINSIZEREL}" )
else ()
  message( FATAL_ERROR "Unknown build type: '${CMAKE_BUILD_TYPE}'" )
endif ()
if ( with-defines )
  foreach ( def ${with-defines} )
    set( ALL_CFLAGS "${def} ${ALL_CFLAGS}" )
    set( ALL_CXXFLAGS "${def} ${ALL_CXXFLAGS}" )
  endforeach ()
endif ()

# all libraries
set( ALL_LIBS
  "-lelements"
  "-lkernel"
  "-llibgrowth"
  "-lio"
  "-lmodels"
  "-lpymodule"
  "-lspatial"
  "${OpenMP_CXX_FLAGS}"
  "${MPI_CXX_LIBRARIES}" )
if ( with-libraries )
  set( ALL_LIBS "${ALL_LIBS};${with-libraries}" )
endif ()
string( REPLACE ";" " " ALL_LIBS "${ALL_LIBS}" )

# all includes
set( ALL_INCLUDES_tmp
  "${PROJECT_SOURCE_DIR}/elements"
  "${PROJECT_SOURCE_DIR}/kernel"
  "${PROJECT_SOURCE_DIR}/libgrowth"
  "${PROJECT_SOURCE_DIR}/io"
  "${PROJECT_SOURCE_DIR}/models"
  "${PROJECT_SOURCE_DIR}/pymodule"
  "${PROJECT_SOURCE_DIR}/spatial"
  "${MPI_CXX_INCLUDE_PATH}" )
set( ALL_INCLUDES "" )
foreach ( INC ${ALL_INCLUDES_tmp} ${with-includes} )
  if ( INC AND NOT INC STREQUAL "" )
    set( ALL_INCLUDES "${ALL_INCLUDES} -I${INC}" )
  endif ()
endforeach ()


###############################################################################
################              Install All Files               #################
###############################################################################

# make a program able to link directly against the library
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
ENDIF()

add_library(cgrowth SHARED ${SOURCE_FILES})
target_link_libraries( cgrowth kernel libgrowth elements spatial )
target_include_directories( cgrowth PUBLIC
    ${PROJECT_SOURCE_DIR}/elements
    ${PROJECT_SOURCE_DIR}/kernel
    ${PROJECT_SOURCE_DIR}/libgrowth
    ${PROJECT_SOURCE_DIR}/io
    ${PROJECT_SOURCE_DIR}/models
    ${PROJECT_SOURCE_DIR}/spatial
)

add_executable(growth "growth.cpp")
target_include_directories( growth PRIVATE
     ${PROJECT_SOURCE_DIR}
     "${PROJECT_SOURCE_DIR}/kernel"
)
target_link_libraries( growth cgrowth )
install(TARGETS growth DESTINATION ${CMAKE_INSTALL_BINDIR})

# this is implemented to prevent python to creat .pyc file
# during developing of pymodule
if (python-develop)
    set(PYTHON_DEVELOPER_MODE "lib64")
endif()

# install target
install(
  TARGETS cgrowth
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# install
install(
  FILES README.md
  DESTINATION "${CMAKE_INSTALL_LIBDIR}${PYEXECDIR}"
)

cgrowth_print_config_summary()
