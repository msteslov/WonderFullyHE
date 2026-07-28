add_library(m2424
    core/abft.cpp
    core/accuracy.cpp
    core/m2424.cpp
    core/profiles.cpp

    ckks/seal_adapter.cpp

    planning/profile_report.cpp
    planning/security_report.cpp
)
add_library(m2424::m2424 ALIAS m2424)

target_include_directories(m2424
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(m2424 PUBLIC cxx_std_17)

target_link_libraries(m2424
    PUBLIC
        SEAL::seal
)
