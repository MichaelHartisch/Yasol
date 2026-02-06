if (NOT HIGHS_ROOT_DIR)
    message(STATUS "No path to the root directory of Highs was provided. I now try to find it myself:")
    file(GLOB_RECURSE HIGHS_DIRS_A "/opt/HiGHS/src/Highs.h")
    file(GLOB_RECURSE HIGHS_DIRS_B "/Applications/HiGHS/src/Highs.h")
    file(GLOB_RECURSE HIGHS_DIRS_C "/opt/HiGHS/highs/Highs.h")
    file(GLOB_RECURSE HIGHS_DIRS_D "/Applications/HiGHS/highs/Highs.h")
    #file(GLOB_RECURSE HiGHS_DIRS_C "**/HiGHS/src/Highs.h")
    list(APPEND HIGHS_DIRS ${HIGHS_DIRS_A})
    list(APPEND HIGHS_DIRS ${HIGHS_DIRS_B})
    list(APPEND HIGHS_DIRS ${HIGHS_DIRS_C})
    list(APPEND HIGHS_DIRS ${HIGHS_DIRS_D})
    if(HIGHS_DIRS)
        message(STATUS "I found these directories:")
        foreach(HIGHS_DIR ${HIGHS_DIRS})
            message(STATUS "${HIGHS_DIR}")
        endforeach()
        foreach(HIGHS_DIR ${HIGHS_DIRS})
            set(HIGHS_ROOT_DIR "${HIGHS_DIR}") #/opt/HiGHS/src/Highs.h
            get_filename_component(HIGHS_ROOT_DIR "${HIGHS_ROOT_DIR}" DIRECTORY)  #/opt/HiGHS/src
            get_filename_component(HIGHS_ROOT_DIR "${HIGHS_ROOT_DIR}" DIRECTORY)  #/opt/HiGHS

            if(NOT EXISTS "${HIGHS_ROOT_DIR}/src/Highs.h" AND NOT EXISTS "${HIGHS_ROOT_DIR}/highs/Highs.h")
                continue()
            else()
                message(STATUS "I did select this folder: ${HIGHS_ROOT_DIR}")
                message(STATUS "If you wish to select another Highs version, specify it. This might look like this: -DHIGHS_ROOT_DIR=/opt/HiGHS")
                break()
            endif()
        endforeach()
    endif()
    if (NOT HIGHS_ROOT_DIR)
        message(FATAL_ERROR "I was not able to find the root directory of Highs. Provide the respective path, which might look like this: -DHIGHS_ROOT_DIR=/opt/HiGHS/src")
   else()
        message(STATUS "Found root. HIGHS_ROOT_DIR: ${HIGHS_ROOT_DIR}")
    endif()
else()
    if(NOT EXISTS "${HIGHS_ROOT_DIR}/src" AND NOT EXISTS "${HIGHS_ROOT_DIR}/highs")
        message(FATAL_ERROR "The provided HIGHS_ROOT_DIR neither seems to contain the folder 'src' nor the folder 'highs'.")
    endif()
endif()

message(STATUS "Looking for HIGHS include files in ${HIGHS_ROOT_DIR}/src and ${HIGHS_ROOT_DIR}/highs and other paths")
find_path(HIGHS_INCLUDE_DIR
  NAMES Highs.h
  HINTS ${HIGHS_ROOT_DIR}/src ${HIGHS_ROOT_DIR}/highs
  PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
  DOC "Path to HiGHS include directory"
)

find_path(HCONFIG_INCLUDE_DIR
  NAMES HConfig.h
  HINTS ${HIGHS_ROOT_DIR}/build
  PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
  DOC "Path to HConfig.h include directory"
)
message(STATUS "HIGHS_INCLUDE_DIR: ${HIGHS_INCLUDE_DIR}")

FIND_LIBRARY(HIGHS_LIBRARY
  NAMES highs
  HINTS ${HIGHS_ROOT_DIR}/build/lib
  PATHS ENV LIBRARY_PATH #unix
        ENV LD_LIBRARY_PATH #unix
  )
message(STATUS "Highs Library: ${HIGHS_LIBRARY}")
