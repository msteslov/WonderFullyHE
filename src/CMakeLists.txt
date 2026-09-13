add_library(m2424
    core/abft.cpp
    core/accuracy.cpp
    core/bootstrap_candidates.cpp
    core/sparse_bootstrap.cpp
    core/bootstrap_contract.cpp
    core/coeff_to_slot_contract.cpp
    core/coeff_to_slot_prefactor.cpp
    core/linear_transform_contract.cpp
    core/m2424.cpp
    core/profiles.cpp

    ckks/evalround_execution.cpp
    ckks/bootstrap.cpp
    ckks/seal_adapter.cpp

    math/canonical_embedding_reference.cpp
    math/coeff_to_slot.cpp
    math/slot_to_coeff.cpp
    math/coeff_to_slot_cert.cpp

    planning/evalround_planner.cpp
    planning/profile_report.cpp
    planning/security_report.cpp
)
add_library(m2424::m2424 ALIAS m2424)

target_include_directories(m2424
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_definitions(m2424 PUBLIC $<$<BOOL:${BUILD_TESTING}>:M2424_ENABLE_SPARSE_DIAGNOSTICS>)
target_compile_features(m2424 PUBLIC cxx_std_17)
target_compile_definitions(m2424 PRIVATE $<$<BOOL:${BUILD_TESTING}>:M2424_ENABLE_MOD_RAISE_CHECKS>)

target_link_libraries(m2424
    PUBLIC
        SEAL::seal
)

option(M2424_ENABLE_EVALMOD_ANALYSIS "Build experimental EvalMod analysis tools" OFF)
if(M2424_ENABLE_EVALMOD_ANALYSIS)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
    pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)
    add_library(m2424_evalmod_analysis
        core/evalmod_domain_analysis.cpp
        math/approximation_lab.cpp
        math/exact_modular_oracle.cpp
        math/evalround_reference.cpp
        math/certified_diagonal.cpp
        planning/slot_to_coeff_preparation.cpp
        planning/coeff_to_slot_preparation.cpp
        planning/sparse_bootstrap.cpp
        planning/bootstrap_planner.cpp
        planning/evalround_execution.cpp
        planning/evalround_synthesis.cpp
        planning/evalmod_cost_model.cpp
        planning/evalmod_feasibility.cpp
        planning/evalmod_synthesis.cpp
    )
    add_library(m2424::evalmod_analysis ALIAS m2424_evalmod_analysis)
    target_include_directories(m2424_evalmod_analysis PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(m2424_evalmod_analysis PUBLIC
        m2424::m2424 PkgConfig::GMPXX PkgConfig::MPFR)
    target_compile_definitions(m2424_evalmod_analysis PRIVATE $<$<BOOL:${BUILD_TESTING}>:M2424_ENABLE_BOOTSTRAP_FIXTURE>)
    target_compile_features(m2424_evalmod_analysis PUBLIC cxx_std_17)
endif()
