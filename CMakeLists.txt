cmake_minimum_required(VERSION 2.8.3)
project(protoHuboGUI)

add_definitions("-Wall -Wextra")

# Modules
set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Eigen
find_package(EIGEN3 3.0.5 QUIET)
if(EIGEN3_FOUND)
  message(STATUS "Looking for EIGEN3 - ${EIGEN3_VERSION} found")
else()
  if(PC_EIGEN3_VERSION)
    message(SEND_ERROR "Looking for EIGEN3 - ${EIGEN3_VERSION} found, ${PROJECT_NAME} requires 3.0.5 or greater.")
  else()
    message(SEND_ERROR "Looking for EIGEN3 - NOT found, please install libeigen3-dev")
  endif()
endif()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

# Boost
set(DART_MIN_BOOST_VERSION 1.46.0 CACHE INTERNAL "Boost min version requirement" FORCE)
if(MSVC OR MSVC90 OR MSVC10)
  add_definitions(-DBOOST_ALL_NO_LIB)
endif()
add_definitions(-DBOOST_TEST_DYN_LINK)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost ${DART_MIN_BOOST_VERSION} COMPONENTS system QUIET)
if(Boost_FOUND)
  message(STATUS "Looking for Boost - ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} found")
else()
  message(SEND_ERROR "Please install system boost version ${DART_MIN_BOOST_VERSION} or higher.")
endif()

# OSG
find_package(OpenSceneGraph 3.0 QUIET
  COMPONENTS osg osgViewer osgManipulator osgGA osgDB)

list(APPEND CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
#list(APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

file(GLOB apps_source "apps/*.cpp")
LIST(SORT apps_source)

message(STATUS "\n    Applications: ")
foreach(app_src_file ${apps_source})
    get_filename_component(app_base ${app_src_file} NAME_WE)
    message(STATUS "Adding application ${app_base}")
    add_executable(${app_base} ${app_src_file})
    target_link_libraries(${app_base} ${exec_dependencies})
endforeach(app_src_file)

target_link_libraries(hubo-ik-gui dart-core dart osgDart ${OPENSCENEGRAPH_LIBRARIES} boost_system assimp HuboCan)
target_link_libraries(hubo-sim dart-core dart osgDart ${OPENSCENEGRAPH_LIBRARIES} boost_system assimp HuboCan)
target_link_libraries(hubo-walk-solver dart-core dart osgDart ${OPENSCENEGRAPH_LIBRARIES} boost_system assimp yaml-cpp HuboCan)
target_link_libraries(hubo-replay dart-core dart osgDart ${OPENSCENEGRAPH_LIBRARIES} boost_system assimp yaml-cpp)
