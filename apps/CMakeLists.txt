add_executable(demo_basic demo_basic.cpp)
target_compile_features(demo_basic PRIVATE cxx_std_17)
target_link_libraries(demo_basic PRIVATE m2424)

set_target_properties(demo_basic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
