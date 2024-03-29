# Check if pymodel changed
set( cppfile "${CMAKE_CURRENT_SOURCE_DIR}/models_manager.cpp" )
set( inputfile "${CMAKE_CURRENT_SOURCE_DIR}/models_manager.cpp.in" )
set( pyfile "${CMAKE_CURRENT_SOURCE_DIR}/models_generator.py" )

if ((NOT EXISTS ${cppfile}) OR (${pyfile} IS_NEWER_THAN ${cppfile}) OR (${inputfile} IS_NEWER_THAN ${cppfile}))
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" "models_generator.py"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
endif()

set( models_sources
    ActinWave.hpp ActinWave.cpp
    models_manager.hpp models_manager.cpp
    # complete model and interfaces
    complete_gc_model.hpp
    extension_interface.hpp
    steering_interface.hpp
    direction_select_interface.hpp
    # extension algorithms
    extension_cst.hpp
    extension_gfluct.hpp extension_gfluct.cpp
    extension_resource_based.hpp extension_resource_based.cpp
    # steering methods
    steering_pull_only.hpp
    steering_memory_based.hpp steering_memory_based.cpp
    steering_srf.hpp steering_srf.cpp
    # direction selection methods
    direction_select_nm.hpp direction_select_nm.cpp
    direction_select_nwa.hpp direction_select_nwa.cpp
    direction_select_rt.hpp direction_select_rt.cpp
    )

add_library( models STATIC ${models_sources} )

set_target_properties(
    models PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME "models_py${Python3_VERSION}"
)

target_include_directories( models PUBLIC
    ${PROJECT_SOURCE_DIR}/elements
    ${PROJECT_SOURCE_DIR}/kernel
    ${PROJECT_SOURCE_DIR}/libgrowth
    ${PROJECT_SOURCE_DIR}/spatial
)
target_link_libraries( models elements kernel libgrowth spatial  )

install( TARGETS models
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
