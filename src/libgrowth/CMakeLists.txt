set( libgrowth_sources
     config.hpp config_impl.hpp config.cpp
     elements_types.hpp
     spatial_types.hpp spatial_types.cpp
     exceptions.hpp exceptions.cpp
     growth_space.hpp growth_space.cpp
     growth_time.hpp growth_time.cpp
     growth_names.hpp growth_names.cpp
)

configure_file(
    ${PROJECT_SOURCE_DIR}/libgrowth/config.hpp.in
    ${PROJECT_SOURCE_DIR}/libgrowth/config.hpp
    @ONLY
)

add_library( libgrowth STATIC ${libgrowth_sources} )

set_target_properties( libgrowth PROPERTIES POSITION_INDEPENDENT_CODE ON )

target_include_directories( libgrowth PRIVATE
    ${PROJECT_SOURCE_DIR}/elements
)

target_link_libraries( libgrowth elements )

install( TARGETS libgrowth
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
