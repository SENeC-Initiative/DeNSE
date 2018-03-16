set( models_sources
    gc_random_walk.hpp  gc_random_walk.cpp
    gc_run_tumble.hpp  gc_run_tumble.cpp
    gc_critical.hpp gc_critical.cpp
    gc_dir_el.hpp gc_dir_el.cpp
    gc_self_referential_forces.hpp gc_self_referential_forces.cpp
    ActinWave.cpp ActinWave.hpp
    models_manager.hpp
    models_manager.cpp
    )

add_library( models STATIC ${models_sources} )

set_target_properties( models PROPERTIES POSITION_INDEPENDENT_CODE ON )

target_include_directories( models PRIVATE
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
