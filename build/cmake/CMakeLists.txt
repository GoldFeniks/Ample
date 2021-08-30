cmake_minimum_required(VERSION 3.12)
project(AMPLE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

if(WIN32)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR};${CMAKE_MODULE_PATH}")

set(PROJECT_DIR ${PROJECT_SOURCE_DIR}/../../)

add_subdirectory(${PROJECT_DIR}/modules/CAMBALA/src/normal_modes ${CMAKE_BINARY_DIR}/normal_modes)
set_target_properties(normal_modes_sample alglib PROPERTIES EXCLUDE_FROM_ALL ON EXCLUDE_FROM_DEFAULT_BUILD ON)

# FFTW3 is currently broken on Linux
if(WIN32)
    find_package(FFTW3 REQUIRED)
    set(FFTW3_LIBRARIES FFTW3::fftw3)
else()
    set(FFTW3_LIBRARIES fftw3)
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

find_package(Threads REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Boost 1.69 COMPONENTS program_options REQUIRED)

include_directories(
    ${PROJECT_DIR}/include 
    ${PROJECT_DIR}/modules/CAMBALA/src/normal_modes
    ${PROJECT_DIR}/modules/delaunay
    ${PROJECT_DIR}/modules/DORK/include
    ${PROJECT_DIR}/modules/zip/include
    ${Boost_INCLUDE_DIR}
    ${PROJECT_DIR}/modules/CAMBALA/third_party/eigen
)

add_executable(AMPLE
        ${PROJECT_DIR}/build/main.cpp
        ${PROJECT_DIR}/include/config.hpp
        ${PROJECT_DIR}/include/initial_conditions.hpp
        ${PROJECT_DIR}/include/io/reader.hpp
        ${PROJECT_DIR}/include/io/writer.hpp
        ${PROJECT_DIR}/include/io/convertors.hpp
        ${PROJECT_DIR}/include/modes.hpp
        ${PROJECT_DIR}/include/rays.hpp
        ${PROJECT_DIR}/include/solver.hpp
        ${PROJECT_DIR}/include/boundary_conditions.hpp
        ${PROJECT_DIR}/include/utils/callback.hpp
        ${PROJECT_DIR}/include/utils/fft.hpp
        ${PROJECT_DIR}/include/utils/interpolation.hpp
        ${PROJECT_DIR}/include/utils/progress_bar.hpp
        ${PROJECT_DIR}/include/utils/types.hpp
        ${PROJECT_DIR}/include/utils/utils.hpp
        ${PROJECT_DIR}/include/utils/verbosity.hpp
        ${PROJECT_DIR}/include/utils/multi_optional.hpp
        ${PROJECT_DIR}/include/utils/join.hpp
        ${PROJECT_DIR}/include/utils/assert.hpp
        ${PROJECT_DIR}/include/utils/dimensions.hpp
)

if(WIN32 AND USE_VCPKG)
    set_target_properties(AMPLE PROPERTIES VS_GLOBAL_VcpkgEnabled true)
endif()

target_link_libraries(AMPLE
    Threads::Threads
    nlohmann_json::nlohmann_json
    ${Boost_LIBRARIES} 
    normal_modes
    ${FFTW3_LIBRARIES}
)
