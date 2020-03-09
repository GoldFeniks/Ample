cmake_minimum_required(VERSION 3.12)
project(ACOUSTIC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR};${CMAKE_MODULE_PATH}")

set(PROJECT_DIR ${PROJECT_SOURCE_DIR}/../../)

find_package(Threads REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Boost 1.69 COMPONENTS program_options REQUIRED)
include_directories(
    ${PROJECT_DIR}/include 
    ${PROJECT_DIR}/modules/Acoustics-at-home/src/cambala
    ${PROJECT_DIR}/modules/delaunay
    ${PROJECT_DIR}/modules/DORK/include
    ${Boost_INCLUDE_DIR}
)

add_executable(ACOUSTIC ${PROJECT_DIR}/build/main.cpp)
target_link_libraries(ACOUSTIC 
    Threads::Threads
    nlohmann_json::nlohmann_json
    ${Boost_LIBRARIES} 
    ${PROJECT_DIR}/modules/Acoustics-at-home/binaries/libnormal_modes.a
)
