#Author: Trent Rabe

cmake_minimum_required(VERSION 3.8.0)

# set default build to debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Release or Debug" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -static-libgcc -static-libstdc++ -pthread -lpthread -Wall" CACHE INTERNAL "")

# SET(CMAKE_EXE_LINKER_FLAGS "$-rpath /apps2/gcc/9.2.0/lib64$")

# get folder name as project name
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})

##############################################
# define cmmake variables for file locations #
##############################################
file(GLOB USER_FILES_1
"./src/driver.cpp"
"./src/PolyEval.cpp"
"./headers/PolyEval.hpp"
"./src/Polynomial.cpp"
"./headers/Polynomial.hpp"
"./src/RandCoeffs.cpp"
"./headers/RandCoeffs.hpp"
"./src/tools.cpp"
"./headers/tools.hpp"
)

# file(GLOB USER_FILES_2
# "./src/driver2.cpp"
# "./src/Functions.cpp"
# "./headers/Functions.hpp"
# "./src/RandChoice.cpp"
# "./headers/RandChoice.hpp"
# "./src/ApplyFuncs.cpp"
# "./headers/ApplyFuncs.hpp"
# "./src/MakeMatrix2.cpp"
# "./headers/MakeMatrix2.hpp"
# "./src/Pixel.cpp"
# "./headers/Pixel.hpp"
# )

# Download and unpack googletest at configure time
# configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
# execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
# 	RESULT_VARIABLE result
# 	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
# if(result)
# 	message(FATAL_ERROR "CMake step for googletest failed: ${result}")
# endif()
# execute_process(COMMAND ${CMAKE_COMMAND} --build .
# 	RESULT_VARIABLE result
# 	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
# if(result)
# 	message(FATAL_ERROR "Build step for googletest failed: ${result}")
# endif()
# # Prevent overriding the parent project's compiler/linker
# # settings on Windows
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# # Add googletest directly to our build. This defines
# # the gtest and gtest_main targets.
# add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
# 	${CMAKE_CURRENT_BINARY_DIR}/googletest-build
# 		EXCLUDE_FROM_ALL)

# create an executable for driver 
add_executable( run_driver ${USER_FILES_1} )
target_link_libraries( run_driver)

# create an executable for driver 
# add_executable( run_driver2 ${USER_FILES_2} )
# target_link_libraries( run_driver2)

# create an executable for DFGF_S1 tests 
# add_executable( run_DFGF_S1_tests ${TEST_FILES_DFGF_S1} ${USER_FILES_2} )
# target_link_libraries( run_DFGF_S1_tests gtest_main)

# create an executable for DFGF_T2 tests 
# add_executable( run_DFGF_T2_tests ${TEST_FILES_DFGF_T2} ${USER_FILES_2} )
# target_link_libraries( run_DFGF_T2_tests gtest_main)

# create an executable for driver tests 
# add_executable( run_driver_tests ${TEST_FILES_DRIVER} ${USER_FILES_1} )
# target_link_libraries( run_driver_tests gtest_main jsoncpp)

# create an executable for RandVec tests 
# add_executable( run_RandVec_tests ${TEST_FILES_RANDVEC} ${RANDVEC_FILES} )
# target_link_libraries( run_RandVec_tests gtest_main)

# create an executable for tools tests 
# add_executable( run_tools_tests ${TEST_FILES_TOOLS} ${TOOLS_FILES} )
# target_link_libraries( run_tools_tests gtest_main)

# create an executable for all tests 
# add_executable( run_all_tests ${TEST_FILES_ALL} ${USER_FILES_1} )
# target_link_libraries( run_all_tests gtest_main jsoncpp)

ENABLE_TESTING()