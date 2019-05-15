# pymodule/CMakeLists.txt

if ( HAVE_PYTHON )
  set( PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYEXECDIR} )
  set( PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYEXECDIR} PARENT_SCOPE)
  set( PYMODULE_DIR ${PROJECT_SOURCE_DIR}/pymodule )

  configure_file(${PYMODULE_DIR}/setup.py.in ${PYMODULE_DIR}/setup.py)
  if (PYTHON_DEVELOPER_MODE)
    install( CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH}:${PYTHON_INSTALL_DIR}:${PYTHON_INSTALL_DIR}/$PY_NAME}\")
    execute_process(
        COMMAND ${PYTHON} ${PYMODULE_DIR}/setup.py
            build -b ${CMAKE_CURRENT_BINARY_DIR}
            develop --prefix=${PYTHON_INSTALL_DIR}
                    --install-scripts=${PYTHON_INSTALL_DIR}
                    --install-data=${PYTHON_INSTALL_DIR}/${PY_NAME}
                    --install-headers=${PYTHON_INSTALL_DIR}
                    --install-lib=${PYTHON_INSTALL_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})"
    )
  else()
    install( CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH}:${PYTHON_INSTALL_DIR}:${PYTHON_INSTALL_DIR}/$PY_NAME}\")
    execute_process(
        COMMAND ${PYTHON} ${PYMODULE_DIR}/setup.py
            build -b ${CMAKE_CURRENT_BINARY_DIR}
            install --prefix=${PYTHON_INSTALL_DIR}
                    --install-scripts=${PYTHON_INSTALL_DIR}
                    --install-data=${PYTHON_INSTALL_DIR}/${PY_NAME}
                    --install-headers=${PYTHON_INSTALL_DIR}
                    --install-lib=${PYTHON_INSTALL_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})"
    )
  endif()
else ()
  message( FATAL_ERROR "Python required")
endif ()
