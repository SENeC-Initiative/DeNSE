set( spatial_sources
     search.hpp search.cpp
)

add_library( spatial STATIC ${spatial_sources} )

set_target_properties( spatial PROPERTIES POSITION_INDEPENDENT_CODE ON )

target_include_directories( spatial PRIVATE
    ${PROJECT_SOURCE_DIR}/libgrowth
    ${PROJECT_SOURCE_DIR}/elements
    ${GEOS_INCLUDE_DIR}
)

target_link_libraries( spatial libgrowth elements ${GEOS_LIBRARY})

install( TARGETS spatial
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)




