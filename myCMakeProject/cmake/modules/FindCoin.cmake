# FindCoinOR.cmake
#
# Usage:
#   cmake -DCOIN_ROOT_DIR=/path/to/coin-or ..
#
# We first try to determine COIN_ROOT_DIR if not given
# and verifies that it exists and contains expected subfolders.

# --------------------------------------------------------------------
# 1. Determine COIN_ROOT_DIR
# --------------------------------------------------------------------
if(NOT COIN_ROOT_DIR)
    # Try environment variables first (optional)
    if(DEFINED ENV{COIN_ROOT_DIR})
        set(COIN_ROOT_DIR "$ENV{COIN_ROOT_DIR}")
        message(STATUS "COIN_ROOT_DIR not set, using environment COIN_ROOT_DIR=${COIN_ROOT_DIR}")
    elseif(DEFINED ENV{COIN_ROOT})
        set(COIN_ROOT_DIR "$ENV{COIN_ROOT}")
        message(STATUS "COIN_ROOT_DIR not set, using environment COIN_ROOT=${COIN_ROOT_DIR}")
    else()
        # Fallback to a hardcoded default
        set(_COIN_DEFAULT_ROOT "/opt/coin-or")
        message(STATUS "COIN_ROOT_DIR not provided — trying default: ${_COIN_DEFAULT_ROOT}")
        set(COIN_ROOT_DIR "${_COIN_DEFAULT_ROOT}")
    endif()
endif()

# Normalize path (handles slashes nicely across platforms)
file(TO_CMAKE_PATH "${COIN_ROOT_DIR}" COIN_ROOT_DIR)

message(STATUS "Using COIN_ROOT_DIR='${COIN_ROOT_DIR}'")

# --------------------------------------------------------------------
# 2. Existence check for COIN_ROOT_DIR
# --------------------------------------------------------------------
if(NOT IS_DIRECTORY "${COIN_ROOT_DIR}")
    message(FATAL_ERROR
        "COIN_ROOT_DIR='${COIN_ROOT_DIR}' does not exist or is not a directory.\n"
        "Please provide a valid path, e.g.:\n"
        "  -DCOIN_ROOT_DIR=/path/to/coin-or"
    )
endif()

# --------------------------------------------------------------------
# 3. Check for expected subdirectories (Clp, Cbc, ...)
# --------------------------------------------------------------------
foreach(_subdir Clp Cbc)
    if(NOT IS_DIRECTORY "${COIN_ROOT_DIR}/${_subdir}")
        message(FATAL_ERROR
            "The provided COIN_ROOT_DIR ('${COIN_ROOT_DIR}') does not contain "
            "the '${_subdir}' subdirectory.\n"
            "Check your COIN_ROOT_DIR setting or coin-or installation."
        )
    endif()
endforeach()



message(STATUS "Looking for COIN include files in ${COIN_ROOT_DIR}")

find_path(COIN_INCLUDE_DIR
  NAMES CbcModel.hpp OsiSolverInterface.hpp OsiClpSolverInterface.hpp
  HINTS ${COIN_ROOT_DIR}/dist/include/coin-or
	${COIN_ROOT_DIR}/dist/include/coin
  	${COIN_ROOT_DIR}/include/coin
        ${COIN_ROOT_DIR}/include/coin-or
  PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
  DOC "Path to Coin include directory"
)
message(STATUS "COIN_INCLUDE_DIR: ${COIN_INCLUDE_DIR}")

# Find each library individually
find_library(COIN_CBC_LIBRARY NAMES Cbc HINTS ${COIN_ROOT_DIR}/dist/lib ${COIN_ROOT_DIR}/lib PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
find_library(COIN_OSICBC_LIBRARY NAMES OsiCbc HINTS ${COIN_ROOT_DIR}/dist/lib ${COIN_ROOT_DIR}/lib PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
find_library(COIN_CGL_LIBRARY NAMES Cgl HINTS ${COIN_ROOT_DIR}/dist/lib ${COIN_ROOT_DIR}/lib PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
find_library(COIN_CLP_LIBRARY NAMES Clp HINTS ${COIN_ROOT_DIR}/dist/lib ${COIN_ROOT_DIR}/lib PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
find_library(COIN_COINUTILS_LIBRARY NAMES CoinUtils HINTS ${COIN_ROOT_DIR}/dist/lib ${COIN_ROOT_DIR}/lib PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
find_library(COIN_OSICLP_LIBRARY NAMES OsiClp HINTS ${COIN_ROOT_DIR}/dist/lib ${COIN_ROOT_DIR}/lib PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
find_library(COIN_OSI_LIBRARY NAMES Osi HINTS ${COIN_ROOT_DIR}/dist/lib ${COIN_ROOT_DIR}/lib PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)

# Display found libraries for debugging
message(STATUS "Coin CBC Library: ${COIN_CBC_LIBRARY}")
message(STATUS "Coin OsiCbc Library: ${COIN_OSICBC_LIBRARY}")
message(STATUS "Coin Cgl Library: ${COIN_CGL_LIBRARY}")
message(STATUS "Coin Clp Library: ${COIN_CLP_LIBRARY}")
message(STATUS "Coin CoinUtils Library: ${COIN_COINUTILS_LIBRARY}")
message(STATUS "Coin OsiClp Library: ${COIN_OSICLP_LIBRARY}")
message(STATUS "Coin Osi Library: ${COIN_OSI_LIBRARY}")
message(STATUS "COIN_ROOT_DIR: ${COIN_ROOT_DIR}")
