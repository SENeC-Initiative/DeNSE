# cmake/ConfigureSummary.cmake
#
# This file is part of DeNSE.
#
# Copyright (C) 2019 The SENeC Initiative
# (adapted from the NEST Initiative scripts)
#
# DENSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# DeNSE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DeNSE.  If not, see <http://www.gnu.org/licenses/>.

function( CGROWTH_PRINT_CONFIG_SUMMARY )
  message( "" )
  message( "--------------------------------------------------------------------------------" )
  message( "DeNSE Configuration Summary" )
  message( "--------------------------------------------------------------------------------" )
  message( "" )
  message( "Build type          : ${CMAKE_BUILD_TYPE}" )
  message( "Target System       : ${CMAKE_SYSTEM_NAME}" )
  message( "Cross Compiling     : ${CMAKE_CROSSCOMPILING}" )
  message( "C compiler          : ${CMAKE_C_COMPILER}" )
  message( "C compiler flags    : ${ALL_CFLAGS}" )
  message( "C++ compiler        : ${CMAKE_CXX_COMPILER}" )
  message( "C++ compiler flags  : ${ALL_CXXFLAGS}" )
  message( "" )
  message( "GEOS library        : Yes ${GEOS_VERSION} ")
  message( "     Includes       : ${GEOS_INCLUDE_DIR} ")
  message( "     Libraries      : ${GEOS_LIBRARY} ")
  message( "BOOST headers       : Yes ${Boost_VERSION} ")
  message( "      Includes      : ${Boost_INCLUDE_DIR} ")
  message( "" )
  message( "Python bindings     : Yes (Python ${PYTHON_VERSION}: ${PYTHON})" )
  message( "       Includes     : ${PYTHON_INCLUDE_DIRS}" )
  message( "       Libraries    : ${PYTHON_LIBRARIES}" )
  message( "" )
  if ( CYTHON_FOUND )
    message( "Cython bindings     : Yes (Cython ${CYTHON_VERSION}: ${CYTHON_EXECUTABLE})" )
  else ()
    message( "Cython bindings     : No. Make sure to cythonize `_pygrowth.pyx` before." )
  endif ()

  if ( OPENMP_FOUND )
    message( "Use threading       : Yes (OpenMP: ${OpenMP_CXX_FLAGS})" )
  else ()
    message( "Use threading       : No" )
  endif ()

  message( "" )

  if ( DOXYGEN_FOUND )
    message( "Use doxygen         : Yes (${DOXYGEN_EXECUTABLE})" )
  else ()
    message( "Use doxygen         : No (no C++ documentation)" )
  endif ()

  if ( SPHINX_FOUND )
    message( "Use Sphinx          : Yes (${SPHINX_EXECUTABLE}) "
             "=> `-Dwith-docs=ON` can be used" )
  else ()
    message( "Use Sphinx          : No (no local documentation)" )
  endif ()

  if ( BREATHE_FOUND )
    message( "Use Breathe         : Yes" )
  else ()
    message( "Use Breathe         : No (no C++ documentation)" )
  endif ()

  message( "" )

  if ( HAVE_MPI )
    message( "Use MPI             : Yes (MPI: ${MPI_CXX_COMPILER})" )
    message( "    FLAGS           : ${MPI_CXX_COMPILE_FLAGS}" )
    message( "    Includes        : ${MPI_CXX_INCLUDE_PATH}" )
    message( "    Link Flags      : ${MPI_CXX_LINK_FLAGS}" )
    message( "    Libraries       : ${MPI_CXX_LIBRARIES}" )
  else ()
    message( "Use MPI             : No" )
  endif ()
  message( "" )

  message( "Python bindings will be installed to:" )
  message( "    ${PYTHON_INSTALL_DIR}" )
  message( "Libraries will be installed to:" )
  message( "    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" )
  message( "Binaries will be installed to:" )
  message( "    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}" )
  message( "" )

  if ( WITH_DOCS )
    message( "Documentation will be installed to ${CMAKE_INSTALL_PREFIX}/${PY_NAME}_doc" )
    message( "" )
  endif ()
  message( "--------------------------------------------------------------------------------" )
  message( "" )

  if ( SUCCESS )
    message( "You can now build and install DeNSE with" )
    if (MSVC)
      message( "  cmake --build . --config Release --target INSTALL" )
      if ( WITH_DOCS )
        message( "  cmake --build . --config Release --target doc (to install documentation)" )
      endif ()
    else ()
      message( "  make" )
      message( "  make install" )
      if ( WITH_DOCS )
        message( "  make doc (to install documentation)" )
      endif ()
    endif()
    message( "" )
  endif ()

endfunction()
