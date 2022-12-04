project(fft_test)
cmake_minimum_required(VERSION 3.24)
set(CXX_STANDARD 11)

include(FetchContent)

add_executable(fft-test fft.cpp)
find_package(Matplot++ REQUIRED)
target_link_libraries(fft-test Matplot++::matplot)
