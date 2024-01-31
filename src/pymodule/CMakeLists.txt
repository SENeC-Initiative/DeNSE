# pymodule/CMakeLists.txt

if ( HAVE_PYTHON )

  set( PYMODULE_DIR ${PROJECT_SOURCE_DIR}/pymodule )

  set( ngexists 0 )
  
  # Check if files changed for developers
  if (python-develop)
    set( cppfile "${PYMODULE_DIR}/dense/_pygrowth.cpp" )
    set( cmakefile "${PYMODULE_DIR}/CMakeLists.txt" )
    file( GLOB_RECURSE files "${PYMODULE_DIR}/*" )
    file( GLOB_RECURSE dirs "${Python3_INSTALL_DIR}/*" )
    set( fileschanged 0 )
    
    foreach(f ${files})
      string(FIND ${f} ".pyx" IS_PYX)
      string(FIND ${f} ".pxd" IS_PXD)
      if ((IS_PYX GREATER -1 OR IS_PXD GREATER -1) AND ${f} IS_NEWER_THAN ${cppfile})
          set( fileschanged 1 )
      endif()
    endforeach(f)
    
    foreach(f ${dirs})
      string(FIND ${f} "dense" NG_INSTALLED)
      if (NOT (NG_INSTALLED EQUAL -1))
          set( ngexists 1 )
      endif()
    endforeach(f)
  endif()
    
  # Compile only if files changed or first install (for developers)
  if ( NOT fileschanged AND ngexists )
      # just move the files (skip the build)
      set( SKIP_BUILD 1 )
  else()
      set( SKIP_BUILD 0 )
  endif()

  # Define MSVC and APPLE variable to 0 for python if they are not defined
  if (NOT MSVC)
    set(MSVC 0)
  endif ()

  if (NOT APPLE)
    set(APPLE 0)
  endif ()

  # configure _pygrowth.pxd and setup.py
  configure_file(${PYMODULE_DIR}/setup.py.in ${PYMODULE_DIR}/setup.py @ONLY)
  configure_file(${PYMODULE_DIR}/dense/_pygrowth.pxd.in ${PYMODULE_DIR}/dense/_pygrowth.pxd @ONLY)

  # run setup.py
  if (MSVC)
    # reformat path for Windows because it is so annoying it does not recognize
    # slash and backslash, also wrap PYTHON between brackets because, hey,
    # there are whitespaces in "Program Files" (of course)!!
    set(STR_CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH};${Python3_INSTALL_DIR}\")")
    install( CODE "${STR_CODE}
      execute_process(
      COMMAND \"${PYTHON}\" ${PYMODULE_DIR}/setup.py
        build -b ${CMAKE_CURRENT_BINARY_DIR}
        install --prefix=${Python3_INSTALL_DIR}
            --install-scripts=${Python3_INSTALL_DIR}
            --install-data=${Python3_INSTALL_DIR}/${PY_NAME}
            --install-headers=${Python3_INSTALL_DIR}
            --install-lib=${Python3_INSTALL_DIR}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})"
    )
  else ()
	  if (Python3_DEVELOPER_MODE)
      install( CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH}:${Python3_INSTALL_DIR}:${Python3_INSTALL_DIR}/${PY_NAME}\")
      execute_process(
        COMMAND ${PYTHON} ${PYMODULE_DIR}/setup.py
          build -b ${CMAKE_CURRENT_BINARY_DIR}
          develop --prefix=${Python3_INSTALL_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})"
      )
	  else()
      install( CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH}:${Python3_INSTALL_DIR}:${Python3_INSTALL_DIR}/${PY_NAME}\")
      execute_process(
        COMMAND ${PYTHON} ${PYMODULE_DIR}/setup.py
          build -b ${CMAKE_CURRENT_BINARY_DIR}
          install --prefix=${Python3_INSTALL_DIR}
              --install-scripts=${Python3_INSTALL_DIR}
              --install-data=${Python3_INSTALL_DIR}/${PY_NAME}
              --install-headers=${Python3_INSTALL_DIR}
              --install-lib=${Python3_INSTALL_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})"
      )
	  endif()
	endif ()
else ()
  message( FATAL_ERROR "Python required")
endif ()
