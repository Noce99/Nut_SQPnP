add_library (SQPnP pnp_solver.cpp)

add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

SET(INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} )

find_package(Eigen3 REQUIRED)
list(APPEND INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} )

include_directories(${INCLUDE_DIRS} )

set(sources test.cpp )
add_executable(test ${sources})

# when other libraries or executables link to SQPnP
target_include_directories (SQPnP PUBLIC ${INCLUDE_DIRS}) 
