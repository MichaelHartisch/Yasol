# FindGurobi.cmake - Locate the Gurobi C API (no C++ wrapper required).
#
# Specify -DGUROBI_ROOT_DIR=/path/to/gurobi1200/linux64 or set GUROBI_HOME.
# A parent directory such as /opt/gurobi1200 is also accepted.
#
# Provides:
#   Gurobi_FOUND, GUROBI_INCLUDE_DIR, GUROBI_LIBRARY, GUROBI_LIBRARY_DIR
#   GUROBI_LIBRARIES, Gurobi::C

include(FindPackageHandleStandardArgs)

set(GUROBI_ROOT_DIR "${GUROBI_ROOT_DIR}" CACHE PATH "Gurobi installation root")

if(NOT GUROBI_ROOT_DIR AND DEFINED ENV{GUROBI_HOME})
    set(GUROBI_ROOT_DIR "$ENV{GUROBI_HOME}")
endif()

set(_GUROBI_PREFIXES)
if(GUROBI_ROOT_DIR)
    list(APPEND _GUROBI_PREFIXES "${GUROBI_ROOT_DIR}")
else()
    # Common Linux and macOS Gurobi installation layouts.
    file(GLOB _GUROBI_INSTALLS LIST_DIRECTORIES true
        "/opt/gurobi*" "/Library/gurobi*" "/usr/local/gurobi*")
    list(SORT _GUROBI_INSTALLS ORDER DESCENDING)
    list(APPEND _GUROBI_PREFIXES ${_GUROBI_INSTALLS})
endif()

set(_GUROBI_INCLUDE_HINTS)
foreach(_prefix IN LISTS _GUROBI_PREFIXES)
    list(APPEND _GUROBI_INCLUDE_HINTS
        "${_prefix}/include"
        "${_prefix}/linux64/include"
        "${_prefix}/macos_universal2/include"
        "${_prefix}/mac64/include"
        "${_prefix}/win64/include")
endforeach()

# An explicitly selected root must not silently fall back to another SDK.
set(_GUROBI_SEARCH_MODE)
if(GUROBI_ROOT_DIR)
    set(_GUROBI_SEARCH_MODE NO_DEFAULT_PATH)
endif()

find_path(GUROBI_INCLUDE_DIR
    NAMES gurobi_c.h
    HINTS ${_GUROBI_INCLUDE_HINTS}
    ${_GUROBI_SEARCH_MODE}
    DOC "Directory containing gurobi_c.h")

# Gurobi's C library suffix encodes major/minor: 12.0 -> gurobi120,
# 13.0 -> gurobi130, etc. Derive it from the header to avoid mixing SDKs.
set(_GUROBI_LIB_NAMES)
if(GUROBI_INCLUDE_DIR)
    file(STRINGS "${GUROBI_INCLUDE_DIR}/gurobi_c.h" _major_line
        REGEX "^[ \t]*#[ \t]*define[ \t]+GRB_VERSION_MAJOR[ \t]+[0-9]+")
    file(STRINGS "${GUROBI_INCLUDE_DIR}/gurobi_c.h" _minor_line
        REGEX "^[ \t]*#[ \t]*define[ \t]+GRB_VERSION_MINOR[ \t]+[0-9]+")
    if(_major_line AND _minor_line)
        string(REGEX REPLACE ".*GRB_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" _major "${_major_line}")
        string(REGEX REPLACE ".*GRB_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" _minor "${_minor_line}")
        math(EXPR _gurobi_lib_version "${_major} * 10 + ${_minor}")
        list(APPEND _GUROBI_LIB_NAMES "gurobi${_gurobi_lib_version}")
    endif()
endif()

# Fall back to known names only if the header version is not readable.
if(NOT _GUROBI_LIB_NAMES)
    list(APPEND _GUROBI_LIB_NAMES gurobi140 gurobi130 gurobi120 gurobi110 gurobi100)
endif()

set(_GUROBI_LIBRARY_HINTS)
# Put the library belonging to the selected header first.
if(GUROBI_INCLUDE_DIR)
    get_filename_component(_gurobi_sdk_root "${GUROBI_INCLUDE_DIR}" DIRECTORY)
    list(APPEND _GUROBI_LIBRARY_HINTS
        "${_gurobi_sdk_root}/lib" "${_gurobi_sdk_root}/lib64")
endif()
foreach(_prefix IN LISTS _GUROBI_PREFIXES)
    list(APPEND _GUROBI_LIBRARY_HINTS
        "${_prefix}/lib" "${_prefix}/lib64"
        "${_prefix}/linux64/lib" "${_prefix}/macos_universal2/lib"
        "${_prefix}/mac64/lib" "${_prefix}/win64/lib")
endforeach()

find_library(GUROBI_LIBRARY
    NAMES ${_GUROBI_LIB_NAMES}
    HINTS ${_GUROBI_LIBRARY_HINTS}
    ${_GUROBI_SEARCH_MODE}
    DOC "Gurobi C shared library (e.g. libgurobi120.so)")

find_package_handle_standard_args(Gurobi
    REQUIRED_VARS GUROBI_INCLUDE_DIR GUROBI_LIBRARY)

if(Gurobi_FOUND)
    get_filename_component(GUROBI_LIBRARY_DIR "${GUROBI_LIBRARY}" DIRECTORY)
    set(GUROBI_LIBRARIES "${GUROBI_LIBRARY}")
    if(NOT TARGET Gurobi::C)
        add_library(Gurobi::C UNKNOWN IMPORTED)
        set_target_properties(Gurobi::C PROPERTIES
            IMPORTED_LOCATION "${GUROBI_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIR}")
    endif()
    message(STATUS "Gurobi headers: ${GUROBI_INCLUDE_DIR}")
    message(STATUS "Gurobi library: ${GUROBI_LIBRARY}")
endif()

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY)
