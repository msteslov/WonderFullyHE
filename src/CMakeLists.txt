add_library(m2424
    core/abft.cpp
    core/accuracy.cpp
    core/bootstrap_candidates.cpp
    core/coeff_to_slot_contract.cpp
    core/evalmod_domain_analysis.cpp
    core/m2424.cpp
    core/profiles.cpp

    ckks/seal_adapter.cpp

    math/canonical_embedding_reference.cpp
    math/coeff_to_slot.cpp
    math/approximation_lab.cpp
    math/exact_modular_oracle.cpp

    planning/profile_report.cpp
    planning/security_report.cpp
    planning/evalmod_cost_model.cpp
)
add_library(m2424::m2424 ALIAS m2424)

find_path(GMP_INCLUDE_DIR gmpxx.h PATHS /opt/homebrew/include /usr/local/include)
find_library(GMP_LIBRARY gmp PATHS /opt/homebrew/lib /usr/local/lib)
find_library(GMPXX_LIBRARY gmpxx PATHS /opt/homebrew/lib /usr/local/lib)
if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP and GMP C++ bindings are required for the exact modular oracle")
endif()

target_include_directories(m2424
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        ${GMP_INCLUDE_DIR}
)

target_compile_features(m2424 PUBLIC cxx_std_17)
target_compile_definitions(m2424 PRIVATE $<$<BOOL:${BUILD_TESTING}>:M2424_ENABLE_MOD_RAISE_CHECKS>)

target_link_libraries(m2424
    PUBLIC
        SEAL::seal
        ${GMPXX_LIBRARY}
        ${GMP_LIBRARY}
)
