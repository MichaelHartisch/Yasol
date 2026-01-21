cmake_minimum_required(VERSION 3.16)
project(myYasol LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")


option(USE_CPLEX "Compile with CPLEX support" OFF)
option(USE_HIGHS "Compile with HiGHS support" OFF)
option(USE_COIN "Compile with COIN-OR CLP/CBC support" OFF)
option(USE_CGL "Compile with COIN-OR CGL to obtain cuts from" FALSE)

set(USE_CPLEX_NUM 0)
set(USE_HIGHS_NUM 0)
set(USE_COIN_NUM 0)

if(USE_CPLEX)
    message(STATUS "Detected architecture: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
    set(USE_CPLEX_NUM 1)
    if(NOT CPLEX_ROOT_DIR)
        message(STATUS "CPLEX_ROOT_DIR is not set. It would be best to specify it with -DCPLEX_ROOT_DIR=/path/to/cplex.")
    else()
        message(STATUS "Using CPLEX_ROOT_DIR: ${CPLEX_ROOT_DIR}")
    endif()
endif()

if(USE_HIGHS)
  set(USE_HIGHS_NUM 1)
endif()

if(USE_COIN)
  set(USE_COIN_NUM 1)
endif()

math(EXPR TOTAL_OPTIONS "${USE_CPLEX_NUM} + ${USE_HIGHS_NUM} + ${USE_COIN_NUM}")

if(TOTAL_OPTIONS GREATER 1)
    message(FATAL_ERROR "More than one solver option is selected. Please select only one.")
elseif(TOTAL_OPTIONS EQUAL 0)
    message(FATAL_ERROR "No solver option is selected. Please select one solver option (USE_CPLEX, USE_HIGHS, or USE_COIN). \n 
    This may look like: cmake -B build -DUSE_CPLEX=ON")
endif()

# Set standard compile options
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -std=c++11 -DIL_STD -m64 -fno-omit-frame-pointer -fno-builtin-malloc -fmessage-length=0")
set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wno-return-type")
set(CMAKE_CXX_FLAGS_DEBUG           "${CMAKE_CXX_FLAGS_DEBUG} -Og -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE         "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3 -flto")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

# Enable sanitizers
option(ENABLE_SANITIZERS "Enable sanitizers" OFF)
if(ENABLE_SANITIZERS)
    message(STATUS "Sanitizers enabled")
    set(SANITIZER_FLAGS "-fsanitize=address" "-fsanitize=undefined")
endif()

# ------------------------------------------------------
# SOURCE FILE COLLECTION (can be replaced with lists)
# ------------------------------------------------------

include_directories(
    ./src
    ${YASOL_PATH}/Yasol/src
    ${YASOL_PATH}/Solver/include
    ${YASOL_PATH}/QBPSolver/src
)


# Add backend-specific sources
if(USE_CPLEX)
    find_package(Cplex REQUIRED)
    include_directories(SYSTEM ${CPLEX_INCLUDE_DIR})
    list(APPEND SOLVER_SRC ${CMAKE_SOURCE_DIR}/Solver/src/ExternSolvers/QpExtSolCplexC.cpp)
    set(EXTERN_SOLVER_INCLUDE_DIR "${CPLEX_INCLUDE_DIR}")
    set(EXTERN_SOLVER_LIB "${CPLEX_LIBRARY}")
endif()

if(USE_HIGHS)
    find_package(Highs REQUIRED)
    include_directories(SYSTEM ${HIGHS_INCLUDE_DIR} ${HCONFIG_INCLUDE_DIR})
    list(APPEND SOLVER_SRC ${CMAKE_SOURCE_DIR}/Solver/src/ExternSolvers/QpExtSolHighs.cpp)
    set(EXTERN_SOLVER_LIB "${HIGHS_LIBRARY}")
    set(EXTERN_SOLVER_INCLUDE_DIR "${HIGHS_INCLUDE_DIR} ${HCONFIG_INCLUDE_DIR}")
endif()

if(USE_COIN)
    find_package(Coin REQUIRED)
    include_directories(SYSTEM ${COIN_INCLUDE_DIR})
    list(APPEND SOLVER_SRC
      ${CMAKE_SOURCE_DIR}/Solver/src/ExternSolvers/QpExtSolCBC.cpp
      ${CMAKE_SOURCE_DIR}/Solver/src/ExternSolvers/QpExtSolCLP.cpp
    )
    set(EXTERN_SOLVER_LIB 
    ${COIN_CBC_LIBRARY}
    ${COIN_OSICBC_LIBRARY}
    ${COIN_CGL_LIBRARY}
    ${COIN_CLP_LIBRARY}
    ${COIN_COINUTILS_LIBRARY}
    ${COIN_OSICLP_LIBRARY}
    ${COIN_OSI_LIBRARY}
    )
    set(EXTERN_SOLVER_INCLUDE_DIR "${COIN_INCLUDE_DIR}")

endif()

if(USE_CGL AND NOT USE_COIN)
    find_package(Coin REQUIRED)
    include_directories(SYSTEM ${COIN_INCLUDE_DIR})
    list(APPEND SOLVER_SRC ${CMAKE_SOURCE_DIR}/Solver/src/ExternSolvers/QpExtSolCBC.cpp)
    set(EXTERN_SOLVER_LIB
    ${EXTERN_SOLVER_LIB}
    ${COIN_CBC_LIBRARY}
    ${COIN_OSICBC_LIBRARY}
    ${COIN_CGL_LIBRARY}
    ${COIN_CLP_LIBRARY} 
    ${COIN_COINUTILS_LIBRARY}
    ${COIN_OSICLP_LIBRARY}
    ${COIN_OSI_LIBRARY}
    )
endif()


# Solver flags
if(USE_CPLEX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCOMPILE_WITH_CPLEX_C -DUSE_CPLEX_C -DUSE_NBD_CPLEX_C")

endif()
if(USE_HIGHS)
    message(STATUS "Will compile with highs.")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HIGHS -DCOMPILE_WITH_HIGHS -DUSE_NBD_HIGHS")

endif()
if(USE_COIN)
    message(STATUS "Will compile with coin-or tools.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CBC -DUSE_NBD_CLP -DCOMPILE_WITH_CLP -DCOMPILE_WITH_CBC")

endif()

add_executable(MyYasol ${CMAKE_SOURCE_DIR}/src/myYasol.cpp)

target_link_libraries(MyYasol 
    PRIVATE
        ${YASOL_PATH}/build/lib/libYasolLib.a
        ${YASOL_PATH}/build/lib/libQBPSolverLib.a
        ${YASOL_PATH}/build/lib/libSolverLib.a
        ${EXTERN_SOLVER_LIB}
        pthread
        stdc++
        dl
)


if (ENABLE_SANITIZERS)
    target_link_options(MyYasol PRIVATE ${SANITIZER_FLAGS})
endif()

