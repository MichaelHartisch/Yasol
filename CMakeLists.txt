cmake_minimum_required(VERSION 3.16)
project(yasol)

enable_testing()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")


option(USE_CPLEX "Compile with CPLEX support" OFF)
option(USE_HIGHS "Compile with HiGHS support" OFF)
option(USE_COIN "Compile with COIN-OR CLP/CBC support" OFF)
option(USE_CGL "Compile with COIN-OR CGL to obtain cuts from" OFF)

set(USE_CPLEX_NUM 0)
set(USE_HIGHS_NUM 0)
set(USE_COIN_NUM 0)

if(USE_CPLEX)
    if(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
        message(STATUS "Detected Apple Silicon (arm64). Forcing build to x86_64.")
        set(CMAKE_OSX_ARCHITECTURES "x86_64")
    else()
        message(STATUS "Detected architecture: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
    endif()
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

# Architecture for CPLEX on macOS ARM
if(USE_CPLEX AND CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
    message(STATUS "Detected Apple Silicon (arm64), forcing x86_64 for CPLEX.")
    set(CMAKE_OSX_ARCHITECTURES "x86_64")
endif()

# Set standard compile options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -std=c++11 -DIL_STD -m64 -DNO_MPI -fno-omit-frame-pointer -fno-builtin-malloc -fmessage-length=0")
set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wno-return-type")
set(CMAKE_CXX_FLAGS_DEBUG           "${CMAKE_CXX_FLAGS_DEBUG} -Og -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE         "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3 -flto")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(LINKER_OPTIONS                  -flto -Wl,--no-as-needed)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build")

if(CMAKE_BUILD_TYPE STREQUAL "debug")
    message(STATUS "Building with custom debug mode (-O1 -g + sanitizers)")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -g")

    # enable sanitizers automatically
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined")
    # set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined")
endif()


if(CMAKE_BUILD_TYPE STREQUAL "debug3")
    message(STATUS "Building with custom debug3 mode (-g3)")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0")
endif()

# Enable sanitizers
option(ENABLE_SANITIZERS "Enable sanitizers" OFF)
if(ENABLE_SANITIZERS)
    message(STATUS "Sanitizers enabled")
    set(SANITIZER_FLAGS
        -fsanitize=address
#        -fsanitize=undefined
        -fno-omit-frame-pointer
        -fno-builtin-malloc
    )
endif()

if(ENABLE_SANITIZERS)
    foreach(tgt SolverLib QBPSolverLib YasolLib yasol)
        target_compile_options(${tgt} PRIVATE ${SANITIZER_FLAGS})
        target_link_options(${tgt} PRIVATE ${SANITIZER_FLAGS})
    endforeach()
endif()
# ------------------------------------------------------
# SOURCE FILE COLLECTION (can be replaced with lists)
# ------------------------------------------------------


file(GLOB_RECURSE SOLVER_SRC 
${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/ExternCode/QLPReader.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/Parser.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/QlpRelaxer.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/QlpStageRelaxer.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/Log.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/QlpConverter.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/QlpSplitter.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Utilities/QlpStageSolver.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/MemoryManagement/MemoryManager.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/numbers/QpDouble.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/numbers/QpRational.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/qlp/Constraint.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/qlp/Qlp.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/qlp/QlpMatCoeff.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/components/QpObjFunc.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/components/QpRhs.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Datastructures/components/QpVar.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Algorithms/qlp2lp/Qlp2Lp.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Algorithms/nd/NbdStageSolver.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver/src/Algorithms/nd/NbdMaster.cpp
)

file(GLOB_RECURSE QBPSOLVER_SRC ${CMAKE_SOURCE_DIR}/QBPSolver/src/*.cpp ${CMAKE_SOURCE_DIR}/QBPSolver/src/**/*.cpp)
file(GLOB_RECURSE YASOL_SRC ${CMAKE_SOURCE_DIR}/Yasol/src/*.cpp ${CMAKE_SOURCE_DIR}/Yasol/src/*.cc)
list(REMOVE_ITEM YASOL_SRC "${CMAKE_SOURCE_DIR}/Yasol/src/Yasol.cpp")

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
    message(STATUS "SEARCH FOR COIN.")
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
    set(EXTERN_SOLVER_INCLUDE_DIR "${COIN_INCLUDE_DIR}")
endif()



# ------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)     # for .a files
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)     # for .so files
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)     # for executables



add_library(SolverLib STATIC ${SOLVER_SRC})
add_library(QBPSolverLib STATIC ${QBPSOLVER_SRC})
add_library(YasolLib STATIC ${YASOL_SRC})

# Includes

target_include_directories(SolverLib PUBLIC ${CMAKE_SOURCE_DIR}/Solver/include ${EXTERN_SOLVER_INCLUDE_DIR})
target_include_directories(QBPSolverLib PUBLIC ${CMAKE_SOURCE_DIR}/QBPSolver/src ${CMAKE_SOURCE_DIR}/Solver/include ${CMAKE_SOURCE_DIR}/Yasol/src)
target_include_directories(YasolLib PUBLIC ${CMAKE_SOURCE_DIR}/Yasol/src ${CMAKE_SOURCE_DIR}/QBPSolver/src ${CMAKE_SOURCE_DIR}/Solver/include ${EXTERN_SOLVER_INCLUDE_DIR}) 

target_link_libraries(SolverLib PRIVATE ${EXTERN_SOLVER_LIB})
target_link_libraries(YasolLib PRIVATE ${EXTERN_SOLVER_LIB})

# Solver flags
if(USE_CPLEX)
    target_compile_definitions(SolverLib PRIVATE -DUSE_CPLEX_C -DCOMPILE_WITH_CPLEX_C -DUSE_NBD_CPLEX_C)
    target_compile_definitions(YasolLib PRIVATE -DUSE_CPLEX_C -DCOMPILE_WITH_CPLEX_C -DUSE_NBD_CPLEX_C)
    set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -DCOMPILE_WITH_CPLEX_C -DUSE_CPLEX_C -DUSE_NBD_CPLEX_C")

endif()
if(USE_HIGHS)
    target_compile_definitions(YasolLib PRIVATE -DUSE_HIGHS -DCOMPILE_WITH_HIGHS -DUSE_NBD_HIGHS)
    target_compile_definitions(SolverLib PRIVATE -DUSE_HIGHS -DCOMPILE_WITH_HIGHS -DUSE_NBD_HIGHS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HIGHS -DCOMPILE_WITH_HIGHS -DUSE_NBD_HIGHS")

endif()
if(USE_COIN)
    message(STATUS "Will compile with coin-or tools.")

    target_compile_definitions(SolverLib PRIVATE -DUSE_CBC -DUSE_NBD_CLP -DCOMPILE_WITH_CLP -DCOMPILE_WITH_CBC)
    target_compile_definitions(YasolLib PRIVATE -DUSE_CBC -DUSE_NBD_CLP -DCOMPILE_WITH_CLP -DCOMPILE_WITH_CBC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CBC -DUSE_NBD_CLP -DCOMPILE_WITH_CLP -DCOMPILE_WITH_CBC")

endif()
if(NOT USE_CGL AND NOT USE_COIN)
    target_compile_definitions(YasolLib PRIVATE -DNO_CGL)
    target_compile_definitions(SolverLib PRIVATE -DNO_CGL)
    target_compile_definitions(QBPSolverLib PRIVATE -DNO_CGL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_CGL")
endif()

# Sanitizers
if(ENABLE_SANITIZERS)
    foreach(lib SolverLib QBPSolverLib YasolLib)
        target_compile_options(${lib} PRIVATE ${SANITIZER_FLAGS})
        target_link_options(${lib} PRIVATE ${SANITIZER_FLAGS})
    endforeach()
endif()

# ------------------------------------------------------
# OPTIONAL EXECUTABLE
# ------------------------------------------------------

add_executable(yasol
    ${CMAKE_SOURCE_DIR}/Yasol/src/Yasol.cpp  # <== ADD YOUR MAIN FILE HERE
)

# target_compile_options(yasol PRIVATE ${CMAKE_CXX_FLAGS})
target_link_libraries(yasol
    PRIVATE
        SolverLib QBPSolverLib YasolLib
        ${EXTERN_SOLVER_LIB}
        pthread
        stdc++
        dl
)


# ------------------------------------------------------
# INSTALL (optional)
# ------------------------------------------------------

# install(TARGETS SolverLib QBPSolverLib YasolLib
#     ARCHIVE DESTINATION lib
#     LIBRARY DESTINATION lib
#     RUNTIME DESTINATION bin
# )
# install(DIRECTORY Solver/include/ DESTINATION include)

file(GLOB EXAMPLE_FILES "${CMAKE_SOURCE_DIR}/examples/*")

foreach(example_file ${EXAMPLE_FILES})
    get_filename_component(example_name "${example_file}" NAME)

    add_test(
        NAME "solver_${example_name}"
        COMMAND $<TARGET_FILE:yasol> "${example_file}"
    )
endforeach()

add_custom_target(check
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS yasol
)

