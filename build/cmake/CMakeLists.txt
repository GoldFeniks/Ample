cmake_minimum_required(VERSION 3.12)
project(ACOUSTIC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR};${CMAKE_MODULE_PATH}")

set(PROJECT_DIR ${PROJECT_SOURCE_DIR}/../../)

# FFTW3 is currently broken on Linux
if(WIN32)
    find_package(FFTW3 REQUIRED)
endif()

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

if(WIN32)
    if (USE_VCPKG)
        set_target_properties(ACOUSTIC PROPERTIES VS_GLOBAL_VcpkgEnabled true)
    endif()
    
    target_link_libraries(ACOUSTIC 
        Threads::Threads
        nlohmann_json::nlohmann_json
        ${Boost_LIBRARIES} 
        ${PROJECT_DIR}/modules/Acoustics-at-home/lib/normal_modes.lib
        FFTW3::fftw3
    )
else()
    target_link_libraries(ACOUSTIC 
        Threads::Threads
        nlohmann_json::nlohmann_json
        ${Boost_LIBRARIES} 
        ${PROJECT_DIR}/modules/Acoustics-at-home/lib/libnormal_modes.a
        fftw3
    )
endif()
