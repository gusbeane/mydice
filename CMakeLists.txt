# Set CMake required version and project name
cmake_minimum_required (VERSION 2.6)
project (DICE)

# Set the version number
set (DICE_VERSION_MAJOR 4)
set (DICE_VERSION_MINOR 8)
set (DICE_VERSION "${DICE_VERSION_MAJOR}.${DICE_VERSION_MINOR}")

# Say hello
MESSAGE(STATUS "Configuring CMake build system for DICE ${DICE_VERSION}")

# Set the master build flag. This flag is used to signal if the build should be
# performed or if it should be aborted.
set (MASTER_BUILD_FLAG true)

# extra path for DICE
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )
MESSAGE( STATUS "CMAKE_MODULE_PATH  = " ${CMAKE_MODULE_PATH} )

# Set up compilation flags through cmake macros
include(SetupFlags)

# contains the full path to the top level directory of your build tree
MESSAGE( STATUS "PROJECT_BINARY_DIR = " ${PROJECT_BINARY_DIR} )
 
# contains the full path to the root of your project source directory,
# i.e. to the nearest directory where CMakeLists.txt contains the PROJECT() command
MESSAGE( STATUS "PROJECT_SOURCE_DIR = " ${PROJECT_SOURCE_DIR} )

MESSAGE( STATUS "CMAKE_BUILD_TYPE = " ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CMAKE_C_FLAGS  = " ${CMAKE_C_FLAGS} )
MESSAGE( STATUS "LIBRARY TYPE     = " ${LIBTYPE} )
 
MESSAGE( STATUS "core OPT     = " ${OPT} )
MESSAGE( STATUS "OSX detected = " ${OSX} )
MESSAGE( STATUS "BUILD_TYPE   = " ${RELEASE} )

MESSAGE( STATUS "C COMPILER     = " ${CMAKE_C_COMPILER})
MESSAGE( STATUS "CXX COMPILER   = " ${CMAKE_CXX_COMPILER})

# Looking for required libraries
FIND_PACKAGE(GSL REQUIRED)
FIND_PACKAGE(FFTW3 REQUIRED)
if (${ENABLE_THREADS} MATCHES "ON")
	FIND_PACKAGE(OpenMP)
	if(OPENMP_FOUND)
		SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	else()
		message(SEND_ERROR "OpenMP not found - skipping building tests")	
	endif(OPENMP_FOUND)
endif (${ENABLE_THREADS} MATCHES "ON")

# ----- Define all build options ----- #

# Shared library option
option (ENABLE_SHARED "Build shared libraries")

# OpenMP threading
option (ENABLE_THREADS "Enable parallel execution" OFF)
mark_as_advanced(ENABLE_THREADS)

# ----- Setup all build options ----- #

# Shared library
if (${ENABLE_SHARED} MATCHES "ON")
	set (BUILD_SHARED_LIBS ON)
	MESSAGE(STATUS "Shared libraries will be built")
endif (${ENABLE_SHARED} MATCHES "ON")

# OpenMP threading, default to no
set (USE_THREADS 0)
if (${ENABLE_THREADS} MATCHES "ON")
	MESSAGE (STATUS "Parallel execution enabled")
	set (USE_THREADS 1)
endif (${ENABLE_THREADS} MATCHES "ON")

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/src/dice_config.h.in ${CMAKE_SOURCE_DIR}/src/dice_config.h )

# Instruct CMake to create shared libraries
#set (BUILD_SHARED_LIBS ON)

# Add the binary directory to the search path so that the DICE
# configuration file will be found
include_directories ("${PROJECT_BINARY_DIR}/src")

# Add the DICE source directory to the build configuration, include files
# first and then the source files
include_directories ("${PROJECT_SOURCE_DIR}/src")
#add_subdirectory (src)


# Setup a header file to hold configuration variables for Starscream
configure_file (
	"${PROJECT_SOURCE_DIR}/src/dice_config.h.in"
	"${PROJECT_SOURCE_DIR}/src/dice_config.h")

# Set the targets
set (LIB dicetools)
set (BIN dice)

# Set the output directories
set (LIBRARY_OUTPUT_PATH lib)
set (EXECUTABLE_OUTPUT_PATH bin)

# ----------------------------------------------------------
# Install SETUP
# ----------------------------------------------------------

mark_as_advanced(CMAKE_INSTALL_PREFIX)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)  # no -DCMAKE_INSTALL_PREFIX=INSTALL_PATH from the command line
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local" CACHE PATH "" FORCE) # default location /home/${USER}/local
ENDIF()

MESSAGE( STATUS "CMAKE_INSTALL_PREFIX =" ${CMAKE_INSTALL_PREFIX})

# Set the default installation directories
set (LIB_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib)
set (INCLUDE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/include)
set (BIN_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/bin)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

if (OSX)
set(CMAKE_MACOSX_RPATH ON) 
endif()
# ----------------------------------------------------------
#
# Set the C files for the library
FILE(GLOB LIB_C_SRC ${PROJECT_SOURCE_DIR}/src/[a-z]*.c)

# Grab the headers
FILE(GLOB LIB_HEADERS src/[a-z]*.h)

# Add the library
add_library (${LIB} ${LIB_C_SRC})

# Prevent "Undefined symbols" error when linking libraries on Darwin
if(OSX)
	set_target_properties(${LIB} PROPERTIES LINK_FLAGS "-undefined suppress -flat_namespace")
endif(OSX)

link_directories(${PROJECT_BINARY_DIR}/lib ${GSL_PATH}/lib ${FFTW3_PATH}/lib)
include_directories(${PROJECT_SOURCE_DIR}/src ${GSL_PATH}/include ${FFTW3_PATH}/include)

# Build the DICE executable
add_executable(${BIN} src/dice.c)

# Link libraries
if(${ENABLE_THREADS} MATCHES "ON")
	target_link_libraries(${BIN} ${LIB} gsl gslcblas fftw3 fftw3_threads m )
else()
	target_link_libraries(${BIN} ${LIB} gsl gslcblas fftw3 m)
endif()


#install (PROGRAMS ${PROJECT_BINARY_DIR}/bin/${BIN} DESTINATION bin)
get_filename_component(exe ${BIN} NAME_WE)
if (OSX)
  set_target_properties(${exe} PROPERTIES INSTALL_RPATH "@loader_path/../lib") 
endif()
install (TARGETS ${exe} RUNTIME DESTINATION bin)

install (FILES ${PROJECT_BINARY_DIR}/lib/lib${LIB}.${LIBEXT} DESTINATION lib)
