
cmake_minimum_required(VERSION 3.18)

project(GPUExactSubgraphCounting LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -O3 -g -Wall -pthread"
)

set(CMAKE_CUDA_FLAGS
    "${CMAKE_CUDA_FLAGS} -O3 -lineinfo"
)

include_directories(
    ${PROJECT_SOURCE_DIR}
)

set(SUBGRAPH_COUNTING_FILES  

    countingalgo.cu
    distributedcountingalgo.cu
)

add_executable(GPUExactSubgraphCounting.out    ${SUBGRAPH_COUNTING_FILES})

find_package(CUDAToolkit REQUIRED)

target_link_libraries(GPUExactSubgraphCounting.out    PUBLIC CUDA::cudart)
