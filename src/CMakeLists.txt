###############################################################################
#################                   Project                   #################
###############################################################################

project ( cgrowth )

set( PY_NAME "dense" )
set( PY_NAME "${PY_NAME}" PARENT_SCOPE )

# Version number.
set (cgrowth_VERSION_MAJOR 0)
set (cgrowth_VERSION_MINOR 1)
set (cgrowth_VERSION_PATCHLEVEL 0)

set (cgrowth_VERSION_MAJOR ${cgrowth_VERSION_MAJOR} PARENT_SCOPE)
set (cgrowth_VERSION_MINOR ${cgrowth_VERSION_MINOR} PARENT_SCOPE)
set( cgrowth_VERSION_PATCHLEVEL ${cgrowth_VERSION_PATCHLEVEL} PARENT_SCOPE)

set( cgrowth_VERSION "${cgrowth_VERSION_MAJOR}.${cgrowth_VERSION_MINOR}.${cgrowth_VERSION_PATCHLEVEL}" )
set( cgrowth_VERSION "${cgrowth_VERSION}" PARENT_SCOPE)

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

set ( ALL_CXX_FLAGS "" PARENT_SCOPE )


###############################################################################
################              Install All Files               #################
###############################################################################

# make a program able to link directly against the library
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
ENDIF()

IF (WIN32 OR CYGWIN)
	SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
	SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON PARENT_SCOPE)
ENDIF()

add_library(cgrowth SHARED ${SOURCE_FILES})
set_target_properties(cgrowth PROPERTIES OUTPUT_NAME "cgrowth_py${PYTHON_VERSION}")
target_link_libraries( cgrowth kernel libgrowth elements spatial ${GEOS_LIBRARY} ${GEOS_C_LIBRARY})
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
     ${PROJECT_SOURCE_DIR}/kernel
     ${PROJECT_SOURCE_DIR}/libgrowth
)

target_link_libraries( growth cgrowth libgrowth )
install(TARGETS growth DESTINATION ${CMAKE_INSTALL_BINDIR})

# install target
install(
  TARGETS cgrowth
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
