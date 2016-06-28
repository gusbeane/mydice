# ============================================================================
# CMake module to detect GSL library
# ============================================================================

SET(GSL_FOUND FALSE)

if ( GSL_PATH ) # user configure cmake with variable -DGSL_PATH="/gsl/path"
	find_library(GSL_LIB NAMES gsl PATHS ${GSL_PATH}/lib NO_DEFAULT_PATH)
	if (GSL_LIB)
		MESSAGE(STATUS "Found GSL_LIB: " ${GSL_LIB})
		SET(GSL_FOUND TRUE)
	endif()
	find_library(GSLCBLAS_LIB NAMES gslcblas PATHS ${GSL_PATH}/lib NO_DEFAULT_PATH)
	if (GSLCBLAS_LIB)
		MESSAGE(STATUS "Found GSLCBLASLIB: " ${GSLCBLAS_LIB})
		SET(GSLCBLAS_FOUND TRUE)
	endif ()
endif ()


if (NOT GSL_FOUND) # try $HOME/local
	find_library(GSL_LIB NAMES gsl PATHS $ENV{HOME}/local/lib)
	find_path(GSL_PATH include/gsl/gsl_math.h HINTS $ENV{HOME}/local)
	if (GSL_PATH)
		MESSAGE(STATUS "Found GSL_PATH: " ${GSL_PATH})
		if(GSL_LIB)
			SET(GSL_FOUND TRUE)
			MESSAGE(STATUS "Found GSL_LIB: " ${GSL_LIB})
		endif()
	endif()
endif()

if (NOT GSLCBLAS_FOUND) # try $HOME/local
	find_library(GSLCBLAS_LIB NAMES gslcblas PATHS $ENV{HOME}/local/lib /opt/local/lib /usr/local/lib)
	if (GSLCBLAS_LIB)
		MESSAGE(STATUS "Found GSLCBLAS_LIB: " ${GSLCBLAS_LIB})
		SET(GSLCBLAS_FOUND TRUE)
	endif()
endif()

if (NOT GSL_FOUND) #  ABORT
	message(SEND_ERROR "GSL not found - skipping building tests")
endif()

if (NOT GSLCBLAS_FOUND) #  ABORT
	message(SEND_ERROR "GSLCBLAS not found - skipping building tests")
endif()

