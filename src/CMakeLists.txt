add_library(m2424
    core/abft.cpp
    core/accuracy.cpp
    core/m2424.cpp
    core/profiles.cpp

    ckks/checked_evaluator.cpp
    ckks/operation_budget.cpp
    ckks/seal_adapter.cpp

    math/diagonal_transform.cpp
    math/linear_transform.cpp
    math/polynomial.cpp

    planning/parameter_planner.cpp
    planning/profile_report.cpp
    planning/security_report.cpp

    research/bootstrap_dft.cpp
    research/bootstrap_layout_v2.cpp
    research/bootstrap_parameters.cpp
    research/bootstrap_pipeline.cpp
    research/bootstrap_plan.cpp
    research/bootstrap_precision_model.cpp
    research/bootstrap_scaling.cpp
    research/bootstrap_stc_reference.cpp
    research/eval_mod.cpp
    research/mod1_approximation.cpp
    research/mod1_circuit.cpp
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
