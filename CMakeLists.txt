set(CMAKE_LEGACY_CYGWIN_WIN32 0)
cmake_minimum_required( VERSION 2.6 )

# add cmake modules: for all `include(...)` first look here
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

if ( COMMAND cmake_policy )
  #  cmake_policy( SET CMP0003 NEW )
  cmake_policy( SET CMP0012 NEW ) # accept "TRUE" "FALSE" in if statements
  cmake_policy( SET CMP0042 NEW ) # MACOSX_RPATH enabled by default
  cmake_policy( SET CMP0054 NEW ) # Only interpret if() arguments as variables or
                                  # keywords when unquoted.
endif()


project ( dense )


set( SUCCESS 1 )


###############################################################################
#################            User-defined options             #################
###############################################################################

set( with-geos ON CACHE STRING "Give directory to GEOS library. [default=ON]" )
set( with-mpi OFF CACHE STRING "Request compilation with MPI. Optionally give directory with MPI installation. [default=OFF]" )
set( with-openmp ON CACHE STRING "Enable OpenMP multithreading. Optional: set OMP flag. [default=ON]" )
set( with-python ON CACHE STRING "Build python bindings. To set a specific Python, set install path. [default=ON]" )
option( cythonize-pybindings "Use Cython to cythonize `_growth.pyx`. If OFF, DeNSE has to be build from a pre-cythonized `_growth.pyx`. [default=ON]" ON )
option( with-docs OFF "Compile documentation (requires `doxygen`, `breathe` and `sphinx`. [default=OFF]" )
option( with-debug OFF "Compile with debug symbols. [default=OFF]" )
option( python-develop OFF "Python compiler allows to modify the code. [default=OFF]" )
set( with-sphinx "" CACHE STRING "Path to sphinx-build. [default=autodetect]" )

# These includes publish function names.
include( ConfigureSummary )
include( GetTriple )
include( ProcessOptions )

# get triples arch-vendor-os
get_host_triple( CGROWTH_HOST_TRIPLE CGROWTH_HOST_ARCH CGROWTH_HOST_VENDOR CGROWTH_HOST_OS )
get_target_triple( CGROWTH_TARGET_TRIPLE CGROWTH_TARGET_ARCH CGROWTH_TARGET_VENDOR CGROWTH_TARGET_OS )

# process compiler and options
if (MSVC)
   set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive- /Zc:twoPhase-" )
else ()
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

    if (CYGWIN AND CMAKE_COMPILER_IS_GNUCC)
       set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U__STRICT_ANSI__" )
    endif ()

    if (CMAKE_BUILD_TYPE MATCHES RELEASE)
       set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-cpp" )
    endif ()
endif ()

cgrowth_process_with_optimize()
cgrowth_process_with_debug()
cgrowth_process_with_warning()
cgrowth_process_with_libraries()
cgrowth_process_with_includes()
cgrowth_process_with_defines()
cgrowth_process_with_python()
cgrowth_process_with_geos()
cgrowth_process_with_boost()
cgrowth_process_with_openmp()
cgrowth_process_with_mpi()
cgrowth_process_with_docs()
cgrowth_get_color_flags()


# check prefix

IF ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
    # set discovered python local folder
    set( PYTHON_INSTALL_DIR "${PY_LOCAL_DIR}" )
    set( PYTHON_INSTALL_DIR "${PY_LOCAL_DIR}" CACHE STRING "python install dir" FORCE)

    # get local user library folder as default install prefix ()
    if (PYTHON_EXECUTABLE MATCHES "Anaconda")
      execute_process(COMMAND conda info --root OUTPUT_VARIABLE CMAKE_INSTALL_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
      string(REPLACE "\\" "/" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
    else ()
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m site --user-base OUTPUT_VARIABLE CMAKE_INSTALL_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif ()

    # sometimes the getuserbase command fails
    if ( "${CMAKE_INSTALL_PREFIX}" STREQUAL "" )
      if (APPLE)
        set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/Library/Python/${PYTHON_VERSION}")
      elseif (MSVC)
        set(CMAKE_INSTALL_PREFIX "\"$ENV{APPDATA}/Python\"")
      else ()
        set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local")
      endif ()
    endif ()

    # cache it
    set( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "install prefix" FORCE)
ELSEIF ( ("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local") AND NOT ("$CACHE{CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local"))
    set( CMAKE_INSTALL_PREFIX "$CACHE{CMAKE_INSTALL_PREFIX}" )
    set( PYTHON_INSTALL_DIR "$CACHE{PYTHON_INSTALL_DIR}" )
ELSEIF ( "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" )
  # /usr/local is forbidden because it is the default CMAKE value and it makes
  # configuration too complex
    message( FATAL_ERROR "Do not install into '/usr/local'. "
             "Please change CMAKE_INSTALL_PREFIX." )
ELSE ()
  # make sure we get the absolute path
  get_filename_component(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)
ENDIF ()

# set libraries in binaries installation directories

if (NOT DEFINED CMAKE_INSTALL_LIBDIR)
  if (UNIX)
    include( GNUInstallDirs )
  elseif (MSVC AND PYTHON_EXECUTABLE MATCHES "Anaconda")
    set(CMAKE_INSTALL_LIBDIR "Lib")
  else ()
    set(CMAKE_INSTALL_LIBDIR "lib")
  endif ()
endif ()

if (NOT DEFINED CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "bin")
endif ()

# check that python install was setup automatically or set it manually

if (NOT DEFINED PYTHON_INSTALL_DIR OR "${PYTHON_INSTALL_DIR}" STREQUAL "")
  set(PYTHON_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${PYEXECDIR}")
  file(MAKE_DIRECTORY "${PYTHON_INSTALL_DIR}")
endif()


###############################################################################
#################        Define Subdirectories here           #################
###############################################################################

add_subdirectory( src )

if ( WITH_DOCS )
  add_subdirectory( docs )
endif ()

configure_file(
    "${PROJECT_SOURCE_DIR}/extra/set_dense_vars.sh.in"
    "${PROJECT_SOURCE_DIR}/set_dense_vars.sh" @ONLY
)


###############################################################################
#################             Summary of flags                #################
###############################################################################

set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CGROWTH_C_COLOR_FLAGS}" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGROWTH_CXX_COLOR_FLAGS}" )

# all compiler flags
if ( NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "None" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS}" )
elseif ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}   ${CMAKE_C_FLAGS_DEBUG}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif ( ${CMAKE_BUILD_TYPE} STREQUAL "Release" )
  set( ALL_CFLAGS "${CMAKE_C_FLAGS}" )
  set( ALL_CXXFLAGS "${CMAKE_CXX_FLAGS}" )
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


###############################################################################
#################               Print summary                 #################
###############################################################################

cgrowth_print_config_summary()
