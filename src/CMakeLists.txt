add_library(m2424
    m2424.cpp
    seal_adapter/seal_adapter.cpp
)

target_include_directories(m2424
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_compile_features(m2424 PUBLIC cxx_std_17)

target_link_libraries(m2424
    PUBLIC
        SEAL::seal
)
