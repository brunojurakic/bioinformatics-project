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
)

add_custom_target(format
    COMMAND clang-format -i ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/src/include/*.h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Formatting source files with clang-format"
)
