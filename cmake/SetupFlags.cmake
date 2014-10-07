if ( NOT SETUP_FLAGS_INCLUDED )

	set( SETUP_FLAGS_INCLUDED 1 )

	macro( typed_cache_set type doc var )
		set ( ${var} ${ARGN} CACHE ${type} ${doc} FORCE )
		set ( ${var} ${ARGN} CACHE ${type} ${doc} FORCE )
	endmacro()

	typed_cache_set ( STRING "setup flags"  SETUP_FLAGS_INCLUDED 1  )

	# Set a default build type if none is given
	if ( NOT CMAKE_BUILD_TYPE ) # Debug default
		typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE Release   )
	else ()                   # else Release
		typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE Debug )
	endif()

	# set optimizer flag
	if ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
		typed_cache_set (STRING "Optimizer" OPT "-g -ggdb")
	else ()
		typed_cache_set (STRING "Optimizer" OPT "-O3")
	endif ()

	typed_cache_set (STRING "compilation warnings" WARNCPP "-Wall")
	typed_cache_set (STRING "compilation warnings" WARNC  "-Wall")
	typed_cache_set (STRING "compilation warnings" WARNF "-Wall" )

	# Just for testing
	option(RELEASE "RELEASE " ON)  # SHARED lib default
	if ( NOT RELEASE ) # Debug TRUE
		typed_cache_set ( BOOL "RELEASE Version" RELEASE FALSE )
	else ()                   # else  FALSE
		typed_cache_set ( BOOL "RELEASE Version" RELEASE TRUE )
	endif()

	# Detect MacOS
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		typed_cache_set ( BOOL "OSX detection" OSX TRUE  )
	else()
		typed_cache_set ( BOOL "OSX detection" OSX FALSE )
	endif()

	# SHARED lib?
	option(BUILD_SHARED_LIBS "Shared lib " ON)  # SHARED lib default
	if ( NOT BUILD_SHARED_LIBS ) # Debug TRUE
		typed_cache_set ( BOOL "Build SHARED" BUILD_SHARED_LIBS FALSE )
	else ()                   # else  FALSE
		typed_cache_set ( BOOL "Build SHARED" BUILD_SHARED_LIBS TRUE)
	endif()

	# Set library type , SHARED or STATIC + LIB extension
	if(BUILD_SHARED_LIBS)
		typed_cache_set ( STRING "LIB TYPE" LIBTYPE SHARED )
		if (OSX)
			typed_cache_set ( STRING "LIB TYPE" LIBEXT  "dylib" )
		else()
			typed_cache_set ( STRING "LIB TYPE" LIBEXT  "so" )
		endif()
	else ()
		typed_cache_set ( STRING "LIB TYPE" LIBTYPE STATIC )
		typed_cache_set ( STRING "LIB TYPE" LIBEXT  "a" )
	endif()

	if (DICE_INSTALLPATH)
		typed_cache_set ( STRING "DICE location" DICEPATH  ${DICE_INSTALLPATH} )
	endif()

endif( NOT SETUP_FLAGS_INCLUDED )
# ============================================================================
