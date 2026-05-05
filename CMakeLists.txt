cmake_minimum_required(VERSION 3.14)
project(gene_variant_finder LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(gene_variant_finder
    src/main.cpp
    src/fastq_parser.cpp
    src/sequence_filter.cpp
    src/fasta_io.cpp
    src/distance.cpp
    src/clustering.cpp
    src/multi_sample.cpp
    src/config.cpp
    src/pipeline.cpp
    src/minimap2_runner.cpp
    src/memory_usage.cpp
)

add_executable(test_runner
    tests/test_main.cpp
    src/distance.cpp
    src/clustering.cpp
    src/sequence_filter.cpp
)
target_include_directories(test_runner PRIVATE ${CMAKE_SOURCE_DIR}/src)

enable_testing()
add_test(NAME unit_tests COMMAND test_runner)

add_custom_target(format
    COMMAND clang-format -i ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/src/include/*.h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Formatting source files with clang-format"
)
