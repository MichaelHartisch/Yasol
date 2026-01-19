
if (NOT CPLEX_ROOT_DIR)
    message(STATUS "No path to the root directory of cplex was provided. I now try to find it myself:")
    file(GLOB_RECURSE CPLEX_DIRS_A "/opt/ibm/ILOG/CPLEX_Studio*/cplex/**/cplex.h")
    file(GLOB_RECURSE CPLEX_DIRS_B "/Applications/CPLEX_Studio*/cplex/**/cplex.h")
    list(APPEND CPLEX_DIRS ${CPLEX_DIRS_A})
    list(APPEND CPLEX_DIRS ${CPLEX_DIRS_B})
    if(CPLEX_DIRS)
        message(STATUS "I found these directories:")
        foreach(CPLEX_DIR ${CPLEX_DIRS})
  	    message(STATUS "${CPLEX_DIR}")
        endforeach()
	foreach(CPLEX_DIR ${CPLEX_DIRS})
            set(CPLEX_ROOT_DIR "${CPLEX_DIR}") #/opt/ibm/ILOG/CPLEX_Studio*/cplex/include/ilcplex/cplex.h
            get_filename_component(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR}" DIRECTORY)  #/opt/ibm/ILOG/CPLEX_Studio*/cplex/include/ilcplex
            get_filename_component(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR}" DIRECTORY)  #/opt/ibm/ILOG/CPLEX_Studio*/cplex/include/ 
            get_filename_component(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR}" DIRECTORY)  #/opt/ibm/ILOG/CPLEX_Studio*/cplex/
            get_filename_component(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR}" DIRECTORY)  #/opt/ibm/ILOG/CPLEX_Studio*

            if(NOT EXISTS "${CPLEX_ROOT_DIR}/cplex")
                continue()
   	    else()
		message(STATUS "I did select this folder: ${CPLEX_ROOT_DIR}")
		message(STATUS "If you wish to select another CPLEX version, specify it. This might look like this: -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio221")
		break()
	    endif()
	endforeach()
    endif()
    if (NOT CPLEX_ROOT_DIR)
        message(FATAL_ERROR "I was not able to find the root directory of cplex. Provide the respective path, which might look like this: -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio*")
#    else()
#        message(STATUS "Found root. CPLEX_ROOT_DIR: ${CPLEX_ROOT_DIR}")
    endif()
else()
    if(NOT EXISTS "${CPLEX_ROOT_DIR}/cplex")
	message(FATAL_ERROR "The provided CPLEX_ROOT_DIR does not seem to contain the folder 'cplex'. Maybe a typo?")
    endif()
endif()
message(STATUS "Looking for CPLEX include files in ${CPLEX_ROOT_DIR}/cplex/include and other paths")

find_path(CPLEX_INCLUDE_DIR
  NAMES ilcplex/cplex.h
  HINTS ${CPLEX_ROOT_DIR}/concert/include
        ${CPLEX_ROOT_DIR}/cplex/include
  PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
  DOC "Path to Cplex include directory"
)
message(STATUS "CPLEX_INCLUDE_DIR: ${CPLEX_INCLUDE_DIR}")

FIND_LIBRARY(CPLEX_LIBRARY
  NAMES cplex
HINTS ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic #linux
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic #osx
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_darwin/static_pic #osx
	${CPLEX_ROOT_DIR}/cplex/lib/arm64_osx/static_pic #osx_arm
  PATHS ENV LIBRARY_PATH #unix
        ENV LD_LIBRARY_PATH #unix
  )
message(STATUS "CPLEX Library: ${CPLEX_LIBRARY}")
if (NOT CPLEX_LIBRARY)
        message(FATAL_ERROR "I was not able to find the libary files of cplex. Provide the respective path which, might look like this: -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio* when calling cmake directly")
else()
	message(STATUS "CPLEX Library: ${CPLEX_LIBRARY}")
endif()
