# pymodule/CMakeLists.txt

if ( HAVE_PYTHON )
  
  set( PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${PYEXECDIR} )
  
  set( PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${PYEXECDIR} PARENT_SCOPE)
  set( PYMODULE_DIR ${PROJECT_SOURCE_DIR}/pymodule )

  # install required python packages
  # execute_process(COMMAND pip ["install --user cython numpy pint scipy shapely"]])

  set( ngexists 0 )
  
  # Check if files changed for developers
  if (python-develop)
    set( cppfile "${PYMODULE_DIR}/dense/_pygrowth.cpp" )
    set( cmakefile "${PYMODULE_DIR}/CMakeLists.txt" )
    file( GLOB_RECURSE files "${PYMODULE_DIR}/*" )
    file( GLOB_RECURSE dirs "${PYTHON_INSTALL_DIR}/*" )
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

  #message( "Build status: ${fileschanged} + ${ngexists}" )
  if (NOT MSVC)
    set(MSVC 0)
  endif ()

  configure_file(${PYMODULE_DIR}/setup.py.in ${PYMODULE_DIR}/setup.py)

  if (MSVC)
    # file(TO_NATIVE_PATH "${PYTHON_INSTALL_DIR}" MS_PYTHON_INSTALL_DIR)
    string(REPLACE "/" "\\\\" MS_PYTHON_INSTALL_DIR "${PYTHON_INSTALL_DIR}")
    set(STR_CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH};${MS_PYTHON_INSTALL_DIR}\")")
	install( CODE "${STR_CODE}
		execute_process(
		COMMAND ${PYTHON} ${PYMODULE_DIR}/setup.py
			build -b ${CMAKE_CURRENT_BINARY_DIR}
			install --prefix=${MS_PYTHON_INSTALL_DIR}
					--install-scripts=${MS_PYTHON_INSTALL_DIR}
					--install-data=${MS_PYTHON_INSTALL_DIR}\\\\${PY_NAME}
					--install-headers=${MS_PYTHON_INSTALL_DIR}
					--install-lib=${MS_PYTHON_INSTALL_DIR}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})"
	)
  else ()
	  if (PYTHON_DEVELOPER_MODE)
		install( CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH}:${PYTHON_INSTALL_DIR}:${PYTHON_INSTALL_DIR}/${PY_NAME}\")
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
		install( CODE "set(ENV{PYTHONPATH} \"$ENV{PYTHONPATH}:${PYTHON_INSTALL_DIR}:${PYTHON_INSTALL_DIR}/${PY_NAME}\")
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
	endif ()
else ()
  message( FATAL_ERROR "Python required")
endif ()
