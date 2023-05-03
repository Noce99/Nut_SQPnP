add_library (SQPnP pnp_solver.cpp)

add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

SET(INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} )

find_package(Eigen3 REQUIRED)
list(APPEND INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} )

find_package(OpenCV)
list(APPEND LIBS ${OpenCV_LIBS})

include_directories(${INCLUDE_DIRS} )

set(sources test.cpp )
add_executable(test ${sources})

set(sources_create_marker create_marker.cpp)
add_executable(create_marker ${sources_create_marker})
target_link_libraries( create_marker ${OpenCV_LIBS} )

set(sources_read_marker read_markers.cpp)
add_executable(read_markers ${sources_read_marker})
target_link_libraries( read_markers ${OpenCV_LIBS} )

set(sources_pnp PnP.cpp)
add_executable(pnp ${sources_pnp})
target_link_libraries( pnp ${OpenCV_LIBS} )

set(sources_pnp PnP_video.cpp)
add_executable(pnp_video ${sources_pnp})
target_link_libraries( pnp_video ${OpenCV_LIBS} )

# when other libraries or executables link to SQPnP
target_include_directories (SQPnP PUBLIC ${INCLUDE_DIRS})   
