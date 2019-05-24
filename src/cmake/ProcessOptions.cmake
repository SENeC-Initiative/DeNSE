# cmake/ProcessOptions.cmake
#
# This file is part of CGROWTH.
#
# Copyright (C) 2019 The SENeC Initiative
# (adapted from the NEST Initiative scripts)
#
# CGROWTH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# CGROWTH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CGROWTH.  If not, see <http://www.gnu.org/licenses/>.

# Here all user defined options will be processed.

# add custom warnings and optimizations
function( CGROWTH_PROCESS_WITH_OPTIMIZE )
  if ( with-optimize )
    if ( with-optimize STREQUAL "ON" )
      set( with-optimize "-O2" )
    else ()
      set( with-optimize "-O1" )
    endif ()
    foreach ( flag ${with-optimize} )
      set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}" PARENT_SCOPE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE )
    endforeach ()
  endif ()
endfunction()

function( CGROWTH_PROCESS_WITH_DEBUG )
  if ( with-debug )
    if ( with-debug STREQUAL "ON" )
      set( with-debug "-g" )
      set( CMAKE_BUILD_TYPE "Debug" PARENT_SCOPE )
    endif ()
    foreach ( flag ${with-debug} )
      set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}" PARENT_SCOPE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE )
    endforeach ()
  else ()
      set( CMAKE_BUILD_TYPE "Release" PARENT_SCOPE )
  endif ()
endfunction()

function( CGROWTH_PROCESS_WITH_WARNING )
  if ( with-warning )
    if ( with-warning STREQUAL "ON" )
      set( with-warning "-Wall" )
    endif ()
    foreach ( flag ${with-warning} )
      set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}" PARENT_SCOPE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE )
    endforeach ()
  endif ()
endfunction()

function( CGROWTH_PROCESS_WITH_LIBRARIES )
  if ( with-libraries )
    if ( with-libraries STREQUAL "ON" )
      message( FATAL_ERROR "-Dwith-libraries requires full library paths." )
    endif ()
    foreach ( lib ${with-libraries} )
      if ( EXISTS "${lib}" )
        link_libraries( "${lib}" )
      else ()
        message( FATAL_ERROR "Library '${lib}' does not exist!" )
      endif ()
    endforeach ()
  endif ()
endfunction()

function( CGROWTH_PROCESS_WITH_INCLUDES )
  if ( with-includes )
    if ( with-includes STREQUAL "ON" )
      message( FATAL_ERROR "-Dwith-includes requires full paths." )
    endif ()
    foreach ( inc ${with-includes} )
      if ( IS_DIRECTORY "${inc}" )
        include_directories( "${inc}" )
      else ()
        message( FATAL_ERROR "Include path '${inc}' does not exist!" )
      endif ()
    endforeach ()
  endif ()
endfunction()

function( CGROWTH_PROCESS_WITH_DEFINES )
  if ( with-defines )
    if ( with-defines STREQUAL "ON" )
      message( FATAL_ERROR "-Dwith-defines requires compiler defines -DXYZ=... ." )
    endif ()
    foreach ( def ${with-defines} )
      if ( "${def}" MATCHES "^-D.*" )
        add_definitions( "${def}" )
      else ()
        message( FATAL_ERROR "Define '${def}' does not match '-D.*' !" )
      endif ()
    endforeach ()
  endif ()
endfunction()

function( CGROWTH_GET_COLOR_FLAGS )
    set( CGROWTH_C_COLOR_FLAGS "" PARENT_SCOPE )
    set( CGROWTH_CXX_COLOR_FLAGS "" PARENT_SCOPE )

    # add colored output from gcc
    if ( CMAKE_C_COMPILER_ID STREQUAL "GNU" )
      if ( NOT CMAKE_C_COMPILER_VERSION VERSION_LESS "4.9" )
        set( CGROWTH_C_COLOR_FLAGS "-fdiagnostics-color=auto" PARENT_SCOPE )
      endif ()
    endif ()
    if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
      if ( NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9" )
        set( CGROWTH_CXX_COLOR_FLAGS "-fdiagnostics-color=auto" PARENT_SCOPE )
      endif ()
    endif ()
endfunction()

function( CGROWTH_PROCESS_WITH_PYTHON )
  # Find Python
  set( HAVE_PYTHON OFF PARENT_SCOPE )
  string(REGEX MATCH "^(2|3)([.][0-9])?" VALID_PYVERSION "${with-python}" )

  if ( ${with-python} STREQUAL "ON" OR VALID_PYVERSION OR EXISTS ${with-python} )

    # Localize the Python interpreter
    if ( ${with-python} STREQUAL "ON" )
      find_package( PythonInterp )
    elseif (EXISTS ${with-python} )
        # Directly get all variables from python
        execute_process(COMMAND "${with-python}" "-c"
          "from distutils import sysconfig as s;import sys;import struct;
print('.'.join(str(v) for v in sys.version_info));
print(sys.prefix);
print(s.get_python_inc(plat_specific=True));
print(s.get_python_lib(plat_specific=True));
print(s.get_config_var('SO'));
print(hasattr(sys, 'gettotalrefcount')+0);
print(struct.calcsize('@P'));
print(s.get_config_var('LDVERSION') or s.get_config_var('VERSION'));
print(s.get_config_var('LIBDIR') or '');
print(s.get_config_var('MULTIARCH') or '');
"
            RESULT_VARIABLE _PYTHON_SUCCESS
            OUTPUT_VARIABLE _PYTHON_VALUES
            ERROR_VARIABLE _PYTHON_ERROR_VALUE)
        # Convert the process output into a list
        if (_PYTHON_SUCCESS MATCHES 0)
          set(PYTHONINTERP_FOUND 1)
        endif ()
        if(WIN32)
          string(REGEX REPLACE "\\\\" "/" _PYTHON_VALUES ${_PYTHON_VALUES})
        endif()
        string(REGEX REPLACE ";" "\\\\;" _PYTHON_VALUES ${_PYTHON_VALUES})
        string(REGEX REPLACE "\n" ";" _PYTHON_VALUES ${_PYTHON_VALUES})
        list(GET _PYTHON_VALUES 0 _PYTHON_VERSION_LIST)
        list(GET _PYTHON_VALUES 1 PYTHON_PREFIX)
        list(GET _PYTHON_VALUES 2 PYTHON_INCLUDE_DIR)
        list(GET _PYTHON_VALUES 3 PYTHON_SITE_PACKAGES)
        list(GET _PYTHON_VALUES 4 PYTHON_MODULE_EXTENSION)
        list(GET _PYTHON_VALUES 5 PYTHON_IS_DEBUG)
        list(GET _PYTHON_VALUES 6 PYTHON_SIZEOF_VOID_P)
        list(GET _PYTHON_VALUES 7 PYTHON_LIBRARY_SUFFIX)
        list(GET _PYTHON_VALUES 8 PYTHON_LIBDIR)
        list(GET _PYTHON_VALUES 9 PYTHON_MULTIARCH)
        # The built-in FindPython didn't always give the version numbers
        string(REGEX REPLACE "\\." ";" _PYTHON_VERSION_LIST ${_PYTHON_VERSION_LIST})
        list(GET _PYTHON_VERSION_LIST 0 PYTHON_VERSION_MAJOR)
        list(GET _PYTHON_VERSION_LIST 1 PYTHON_VERSION_MINOR)
        list(GET _PYTHON_VERSION_LIST 2 PYTHON_VERSION_PATCH)

        # Make sure all directory separators are '/'
        string(REGEX REPLACE "\\\\" "/" PYTHON_PREFIX ${PYTHON_PREFIX})
        string(REGEX REPLACE "\\\\" "/" PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR})
        string(REGEX REPLACE "\\\\" "/" PYTHON_SITE_PACKAGES ${PYTHON_SITE_PACKAGES})
    else ()
      find_package( PythonInterp ${with-python} REQUIRED )
    endif ()

    if ( PYTHONINTERP_FOUND )
      set( PYTHONINTERP_FOUND "${PYTHONINTERP_FOUND}" PARENT_SCOPE )
      set( PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} PARENT_SCOPE )
      set( PYTHON ${PYTHON_EXECUTABLE} PARENT_SCOPE )
      set( PYTHON_VERSION_MAJOR ${PYTHON_VERSION_MAJOR} PARENT_SCOPE )
      set( PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} PARENT_SCOPE )

      # Localize Python lib/header files and make sure that their version matches
      # the Python interpreter version
      if (${PYTHON_VERSION})
        find_package( PythonLibs ${PYTHON_VERSION} EXACT )
      else ()
        find_package( PythonLibs )
      endif ()

      if ( PYTHONLIBS_FOUND )
        set( HAVE_PYTHON ON PARENT_SCOPE )
        # export found variables to parent scope
        set( PYTHONLIBS_FOUND "${PYTHONLIBS_FOUND}" PARENT_SCOPE )
        set( PYTHON_INCLUDE_DIRS "${PYTHON_INCLUDE_DIRS}" PARENT_SCOPE )
        set( PYTHON_LIBRARIES "${PYTHON_LIBRARIES}" PARENT_SCOPE )

        if ( cythonize-pybindings )
          find_package( Cython )
          if ( CYTHON_FOUND )
            # confirmed working: 0.19.2+
            if ( CYTHON_VERSION VERSION_LESS "0.19.2" )
              message( FATAL_ERROR "Your Cython version is too old. Please install "
                                   "newer version (0.19.2+)" )
              set( SUCCESS 0 PARENT_SCOPE )
            endif ()

            # export found variables to parent scope
            set( CYTHON_FOUND "${CYTHON_FOUND}" PARENT_SCOPE )
            set( CYTHON_EXECUTABLE "${CYTHON_EXECUTABLE}" PARENT_SCOPE )
            set( CYTHON_VERSION "${CYTHON_VERSION}" PARENT_SCOPE )
          endif ()
        endif ()

        # set execdir for python files
        if (MSVC AND PYTHON_EXECUTABLE MATCHES "Anaconda")
          set( PYEXECDIR "site-packages" PARENT_SCOPE )
        else ()
          set( PYEXECDIR "python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages" PARENT_SCOPE )
        endif ()
      else ()
        message(
          FATAL_ERROR "Python libraries not found, you requested "
          "${PYTHON_VERSION_STRING} but the associated library cannot be found."
        )
      endif ()
    else ()
        message(
          FATAL_ERROR "Python executable not found (you requested "
          "Python ${with-python})."
        )
    endif ()
  else ()
    message(
      FATAL_ERROR "Invalid option: -Dwith-python=" ${with-python}
      ". Valid options are -Dwith-python=ON, -Dwith-python=2, "
      "-Dwith-python=3, or -Dwith-python=X.Y to chose a specific version."
    )
    set( SUCCESS 0 PARENT_SCOPE )
  endif ()
endfunction()


function( CGROWTH_PROCESS_WITH_GEOS )
  # Find GEOS
  if ( NOT "${with-geos}" STREQUAL "ON" )
      find_package( GEOS REQUIRED PATH ${with-geos} )
  else()
      find_package( GEOS REQUIRED )
  endif ()
  if ( GEOS_FOUND )
    # export found to parent
    set ( GEOS_FOUND ${GEOS_FOUND} PARENT_SCOPE )
  else ()
    message( FATAL_ERROR "GEOS required" )
  endif ()
endfunction()


function( CGROWTH_PROCESS_WITH_BOOST )
  # Find BOOST
  find_package(Boost 1.62 REQUIRED)
  #~ find_package(Boost REQUIRED)
  if ( Boost_FOUND )
    # export found to parent
    set ( Boost_FOUND ${Boost_FOUND} PARENT_SCOPE )
    set ( Boost_VERSION "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}" )
    set ( Boost_VERSION "${Boost_VERSION}" PARENT_SCOPE )

    # use boost hash in functional (< 1.67) or in container_hash (>= 1.67)
    if (Boost_VERSION VERSION_GREATER_EQUAL 1.67)
      set(BOOST_1_67_PLUS ON PARENT_SCOPE)
      # add_compile_definitions(BOOST_1_67_PLUS)
    endif ()
  else ()
    message( FATAL_ERROR "BOOST required" )
  endif ()
endfunction()


function( CGROWTH_PROCESS_WITH_OPENMP )
  # Find OPENMP
  if ( with-openmp )
    if ( NOT "${with-openmp}" STREQUAL "ON" )
      message( STATUS "Set OpenMP argument: ${with-openmp}")
      # set variablesin this scope
      set( OPENMP_FOUND ON )
      set( OpenMP_C_FLAGS "${with-openmp}" )
      set( OpenMP_CXX_FLAGS "${with-openmp}" )
    else ()
      find_package( OpenMP )
    endif ()
    if ( OPENMP_FOUND )
      # export found variables to parent scope
      set( WITH_OMP ON PARENT_SCOPE )
      set( OPENMP_FOUND "${OPENMP_FOUND}" PARENT_SCOPE )
      set( OpenMP_C_FLAGS "${OpenMP_C_FLAGS}" PARENT_SCOPE )
      set( OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS}" PARENT_SCOPE )
      # set flags
      set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" PARENT_SCOPE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" PARENT_SCOPE )
    endif ()
  endif ()
endfunction()

function( CGROWTH_PROCESS_WITH_MPI )
  # Find MPI
  set( WITH_MPI OFF PARENT_SCOPE )
  if ( with-mpi )
    find_package( MPI )
    if ( MPI_CXX_FOUND )
      set( WITH_MPI ON PARENT_SCOPE )

      set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}" PARENT_SCOPE )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}" PARENT_SCOPE )

      set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}" PARENT_SCOPE )
      include_directories( ${MPI_CXX_INCLUDE_PATH} )
      # is linked in nestkernel/CMakeLists.txt

      # export found variables to parent scope
      set( MPI_C_FOUND "${MPI_C_FOUND}" PARENT_SCOPE )
      set( MPI_C_COMPILER "${MPI_C_COMPILER}" PARENT_SCOPE )
      set( MPI_C_COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}" PARENT_SCOPE )
      set( MPI_C_INCLUDE_PATH "${MPI_C_INCLUDE_PATH}" PARENT_SCOPE )
      set( MPI_C_LINK_FLAGS "${MPI_C_LINK_FLAGS}" PARENT_SCOPE )
      set( MPI_C_LIBRARIES "${MPI_C_LIBRARIES}" PARENT_SCOPE )
      set( MPI_CXX_FOUND "${MPI_CXX_FOUND}" PARENT_SCOPE )
      set( MPI_CXX_COMPILER "${MPI_CXX_COMPILER}" PARENT_SCOPE )
      set( MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}" PARENT_SCOPE )
      set( MPI_CXX_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}" PARENT_SCOPE )
      set( MPI_CXX_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}" PARENT_SCOPE )
      set( MPI_CXX_LIBRARIES "${MPI_CXX_LIBRARIES}" PARENT_SCOPE )
      set( MPIEXEC "${MPIEXEC}" PARENT_SCOPE )
      set( MPIEXEC_NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}" PARENT_SCOPE )
      set( MPIEXEC_PREFLAGS "${MPIEXEC_PREFLAGS}" PARENT_SCOPE )
      set( MPIEXEC_POSTFLAGS "${MPIEXEC_POSTFLAGS}" PARENT_SCOPE )
    endif ()
  endif ()
endfunction()

function( FIND_PYTHON_MODULE module )
  string( TOUPPER ${module} module_upper )
  execute_process( COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import ${module}"
    RESULT_VARIABLE _module_status
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if ( NOT _module_status )
    set( PY_${module_upper}_FOUND ON PARENT_SCOPE )
    set( PY_${module_upper} "Yes" )
  else ()
    set( PY_${module_upper}_FOUND OFF PARENT_SCOPE)
    set( PY_${module_upper} "No" )
  endif ()
  find_package_handle_standard_args(
    ${module} DEFAULT_MSG PY_${module_upper}
  )
endfunction()

function( CGROWTH_PROCESS_WITH_DOCS )
  # Find documentation software
  set( WITH_DOCS OFF PARENT_SCOPE )

  find_package( Doxygen )
  if ( DOXYGEN_FOUND )
    set( DOXYGEN_FOUND ON PARENT_SCOPE )
    set( DOXYGEN_EXECUTABLE "${DOXYGEN_EXECUTABLE}" PARENT_SCOPE )
  endif ()

  set( BREATHE_FOUND OFF PARENT_SCOPE )
  set( BREATHE_FOUND OFF )
  if ( HAVE_PYTHON )
    find_python_module(breathe)
    if ( PY_BREATHE_FOUND )
      set( BREATHE_FOUND ON PARENT_SCOPE )
      set( BREATHE_FOUND ON )
    endif ()
  endif ()

  if (BREATHE_FOUND AND DOXYGEN_FOUND)
    set( CPP_DOCS_AVAILABLE True PARENT_SCOPE)
  else ()
    set( CPP_DOCS_AVAILABLE False PARENT_SCOPE)
  endif ()

  if ( "${with-sphinx}" STREQUAL "" )
    find_package( Sphinx )
  else ()
    set( SPHINX_FOUND ON )
    set( SPHINX_FOUND ON PARENT_SCOPE )
    get_filename_component(SPHINX_EXECUTABLE "${with-sphinx}" ABSOLUTE)
  endif ()

  if ( SPHINX_FOUND )
    set( SPHINX_FOUND ON PARENT_SCOPE )
    set( SPHINX_EXECUTABLE "${SPHINX_EXECUTABLE}" PARENT_SCOPE )

    if ( NOT DEFINED SPHINX_THEME )
      set(SPHINX_THEME default)
    endif ()

    if ( NOT DEFINED SPHINX_THEME_DIR )
      set( SPHINX_THEME_DIR )
    endif()
  endif ()

  if ( with-docs AND SPHINX_FOUND )
      set ( WITH_DOCS ON PARENT_SCOPE )
  endif ()
endfunction()
