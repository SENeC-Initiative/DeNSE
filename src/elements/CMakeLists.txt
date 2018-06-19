set( elements_sources
     Branch.hpp Branch.cpp
     Branching.hpp Branching.cpp
     Neurite.hpp Neurite.cpp
     Neuron.hpp Neuron.cpp
     GrowthCone.hpp  GrowthCone.cpp
     Node.hpp       Node.cpp
)

add_library( elements STATIC ${elements_sources} )

set_target_properties( elements PROPERTIES POSITION_INDEPENDENT_CODE ON )

target_include_directories(
    elements PUBLIC
    ${PROJECT_SOURCE_DIR}/kernel
    ${PROJECT_SOURCE_DIR}/libgrowth
    ${PROJECT_SOURCE_DIR}/spatial
    ${PROJECT_SOURCE_DIR}/models
    ${PROJECT_SOURCE_DIR}/io
)

target_link_libraries( elements io models kernel libgrowth spatial )

install( TARGETS elements
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
