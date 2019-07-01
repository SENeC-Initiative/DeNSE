set( kernel_sources
     kernel_manager.hpp kernel_manager.cpp
     manager_interface.hpp
     neuron_manager.hpp neuron_manager.cpp
     parallelism_manager.hpp parallelism_manager.cpp
     recorders.hpp recorders.cpp
     record_manager.hpp record_manager.cpp
     rng_manager.hpp rng_manager.cpp
     simulation_manager.hpp simulation_manager.cpp
     space_manager.hpp space_manager.cpp
)

set( WITH_OMP ${OPENMP_FOUND} )

add_library( kernel STATIC ${kernel_sources} )

target_include_directories(
    kernel PRIVATE
    ${PROJECT_SOURCE_DIR}/elements
    ${PROJECT_SOURCE_DIR}/libgrowth
    ${PROJECT_SOURCE_DIR}/models
    ${PROJECT_SOURCE_DIR}/spatial
    ${PROJECT_SOURCE_DIR}/io
    ${GEOS_INCLUDE_DIR}
)

set_target_properties(
    kernel PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME "kernel_py${PYTHON_VERSION}"
)

target_link_libraries( kernel elements libgrowth models spatial io ${GEOS_LIBRARY})

install( TARGETS kernel
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
