add_library(m2424
    core/abft.cpp
    core/accuracy.cpp
    core/bootstrap_candidates.cpp
    core/coeff_to_slot_contract.cpp
    core/m2424.cpp
    core/profiles.cpp

    ckks/seal_adapter.cpp

    math/canonical_embedding_reference.cpp
    math/coeff_to_slot.cpp

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
target_compile_definitions(m2424 PRIVATE $<$<BOOL:${BUILD_TESTING}>:M2424_ENABLE_MOD_RAISE_CHECKS>)

target_link_libraries(m2424
    PUBLIC
        SEAL::seal
)
