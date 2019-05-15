# docs/CMakeLists.txt

if ( WITH_DOCS )
  # configured documentation tools and intermediate build results
  set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

  # HTML output directory
  set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

  if ( NOT SPHINX_HTML_THEME )
    set ( SPHINX_HTML_THEME "bizstyle" )
    set ( WITH_BOOTSTRAP "True" )
  else ( )
    set ( WITH_BOOTSTRAP "False" )
  endif ()

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
    @ONLY
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
    @ONLY
  )

  add_custom_target( doc
    COMMAND ${CMAKE_COMMAND} -E
        make_directory $BINARY_BUILD_DIR
    COMMAND ${SPHINX_EXECUTABLE}
        -q -b html
        -c ${CMAKE_CURRENT_BINARY_DIR}
        -d ${SPHINX_CACHE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${SPHINX_HTML_DIR}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Compiling documentation..."
  )

  install(CODE "execute_process(COMMAND
    ${SPHINX_EXECUTABLE}
    -q -b html
    -c ${CMAKE_CURRENT_BINARY_DIR}
    -d ${SPHINX_CACHE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SPHINX_HTML_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})"
  )

  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html/" DESTINATION "${CMAKE_INSTALL_PREFIX}/${PY_NAME}_doc")

endif ()
