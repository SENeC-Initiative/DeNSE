set( io_sources
    Skeleton.cpp Skeleton.hpp
    Swc.cpp    Swc.hpp
    )

add_library( io STATIC ${io_sources} )

set_target_properties( io PROPERTIES POSITION_INDEPENDENT_CODE ON )

target_include_directories( io PRIVATE
    ${PROJECT_SOURCE_DIR}/elements
    ${PROJECT_SOURCE_DIR}/libgrowth
    ${PROJECT_SOURCE_DIR}/spatial


)
target_link_libraries( io elements libgrowth spatial  )

install( TARGETS io
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
