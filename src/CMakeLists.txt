add_library(m2424
    accuracy.cpp
    abft.cpp
    bootstrap.cpp
    bootstrap_plan.cpp
    bootstrap_prototype.cpp
    bootstrap_scaling.cpp
    checked_evaluator.cpp
    diagonal_transform.cpp
    eval_mod.cpp
    linear_transform.cpp
    m2424.cpp
    polynomial.cpp
    profiles.cpp
    profile_report.cpp
    security_report.cpp
    seal_adapter/seal_adapter.cpp
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
