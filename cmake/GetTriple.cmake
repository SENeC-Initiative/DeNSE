# Based on libcxx/cmake/Modules/GetTriple.cmake

# Define functions to get the host and target triple.

function( get_host_triple out out_arch out_vendor out_os )
  # Get the architecture.
  set( arch "${CMAKE_HOST_SYSTEM_PROCESSOR}" )
  # i686 is an enhanced version of x86
  if ( arch STREQUAL "x86" )
    set( arch "i686" )
  endif ()

  # Get the vendor.
  if ( "${CMAKE_HOST_SYSTEM_NAME}" MATCHES "^Darwin.*" )
    set( vendor "apple" )
  else ()
    set( vendor "pc" )
  endif ()

  # Get os.
  if ( "${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Windows" )
    set( os "win32" )
  else ()
    string( TOLOWER "${CMAKE_HOST_SYSTEM_NAME}" os )
  endif ()

  set( triple "${arch}-${vendor}-${os}" )

  # return values
  set( ${out} ${triple} PARENT_SCOPE )
  set( ${out_arch} ${arch} PARENT_SCOPE )
  set( ${out_vendor} ${vendor} PARENT_SCOPE )
  set( ${out_os} ${os} PARENT_SCOPE )

  message( STATUS "Host triple: ${triple}" )
endfunction ()


function( get_target_triple out out_arch out_vendor out_os )
  # Get the architecture.
  set( arch "${CMAKE_SYSTEM_PROCESSOR}" )
  # i686 is an enhanced version of x86
  if ( arch STREQUAL "x86" )
    set( arch "i686" )
  endif ()

  # Get the vendor.
  if ( "${CMAKE_SYSTEM_NAME}" MATCHES "^Darwin.*")
    set( vendor "apple" )
  elseif ( TRIPLET_VENDOR )
    # In our own tool chain files we define vendors.
    set( vendor "${TRIPLET_VENDOR}" )
  else ()
    set( vendor "pc" )
  endif ()

  # Get os.
  if ( "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows" )
    set( os "win32" )
  else ()
    string( TOLOWER "${CMAKE_SYSTEM_NAME}" os )
  endif ()

  set(triple "${arch}-${vendor}-${os}")

  # return values
  set( ${out} ${triple} PARENT_SCOPE )
  set( ${out_arch} ${arch} PARENT_SCOPE )
  set( ${out_vendor} ${vendor} PARENT_SCOPE )
  set( ${out_os} ${os} PARENT_SCOPE )

  message(STATUS "Target triple: ${triple}")
endfunction()
