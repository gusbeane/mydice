# ============================================================================
# CMake module to detect FFTW3 library
# ============================================================================

SET(FFTW3_FOUND FALSE)
if(FFTW3_PATH) # user configure cmake with variable -DFFTW3_PATH="/fftw3/path"
	find_library(FFTW3_LIB NAMES fftw3 PATHS ${FFTW3_PATH}/lib NO_DEFAULT_PATH)
	if (FFTW3_LIB)
		MESSAGE(STATUS "Found FFTW3_LIB: " ${FFTW3_LIB})
		SET(FFTW3_FOUND TRUE)
	endif()
endif ()
if(NOT FFTW3_FOUND) # try some standard installations
	find_library(FFTW3_LIB NAMES fftw3 PATHS $ENV{HOME}/local/lib)
	find_path(FFTW3_PATH include/fftw3.h HINTS $ENV{HOME}/local)
	if (FFTW3_PATH)
		MESSAGE(STATUS "Found FFTW3_PATH: " ${FFTW3_PATH})
		if(FFTW3_LIB)
			SET(FFTW3_FOUND TRUE)
			MESSAGE(STATUS "Found FFTW3_LIB: " ${FFTW3_LIB})
		endif()
	endif()
endif()
if (NOT FFTW3_FOUND) #  ABORT
	message(SEND_ERROR "FFTW3 not found - skipping building tests")
endif()

# If threading is activated
if(ENABLE_THREADS MATCHES "ON")
	SET(FFTW3_THREADS_FOUND FALSE)
        if(FFTW3_THREADS_PATH) # user configure cmake with variable -DFFTW3_PATH="/fftw3/path"
                find_library(FFTW3_THREADS_LIB NAMES fftw3_threads PATHS ${FFTW3_THREADS_PATH}/lib NO_DEFAULT_PATH)
                if (FFTW3_THREADS_LIB)
                        MESSAGE(STATUS "Found FFTW3_THREADS_LIB: " ${FFTW3_THREADS_LIB})
                        SET(FFTW3_THREADS_FOUND TRUE)
                endif()
        endif()
        if(NOT FFTW3_THREADS_FOUND) # try $HOME/local
                find_library(FFTW3_THREADS_LIB NAMES fftw3_threads PATHS $ENV{HOME}/local/lib)
		
                if (FFTW3_THREADS_LIB)
                        MESSAGE(STATUS "Found FFTW3_THREADS_LIB: " ${FFTW3_THREADS_LIB})
                        SET(FFTW3_THREADS_FOUND TRUE)
                endif()
        endif()
        if (NOT FFTW3_THREADS_FOUND) #  ABORT
                message(SEND_ERROR "FFTW3_THREADS not found - skipping building tests")
        endif()
endif() 

