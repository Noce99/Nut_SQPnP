add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

SET(INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} )

find_package(Eigen3 REQUIRED)
list(APPEND INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} )

find_package(OpenCV)
list(APPEND LIBS ${OpenCV_LIBS})

include_directories(${INCLUDE_DIRS} )

set(sources_create_marker create_marker.cpp)
add_executable(create_marker ${sources_create_marker})
target_link_libraries( create_marker ${OpenCV_LIBS} )

set(sources_pnp PnP_photo.cpp)
add_executable(pnp_photo ${sources_pnp})
target_link_libraries( pnp_photo ${OpenCV_LIBS} )

set(sources_pnp PnP_video.cpp)
add_executable(pnp_video ${sources_pnp})
target_link_libraries( pnp_video ${OpenCV_LIBS} )
