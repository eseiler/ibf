cmake_minimum_required (VERSION 3.8)

add_executable (build_ibf build_ibf.cpp)
target_link_libraries (build_ibf PUBLIC seqan3::seqan3)

add_executable (query_ibf query_ibf.cpp)
target_link_libraries (query_ibf PUBLIC seqan3::seqan3)
