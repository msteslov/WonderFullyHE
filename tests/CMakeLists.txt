add_library(m2424_test_support STATIC seal_adapter_test_support.cpp)
target_compile_features(m2424_test_support PRIVATE cxx_std_17)
target_compile_definitions(m2424_test_support PRIVATE
    M2424_ENABLE_MOD_RAISE_CHECKS
    M2424_ENABLE_SPARSE_DIAGNOSTICS
)
target_link_libraries(m2424_test_support PUBLIC m2424)

function(m2424_add_test target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424_test_support)
    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES TIMEOUT 60)
endfunction()

set(M2424_ADAPTER_AND_SECURITY_TESTS
    test_accuracy
    test_adapter_failures
    test_bootstrap_candidates
    test_bootstrap_contract
    test_ckks_primitives
    test_sparse_security
)

set(M2424_LINEAR_TRANSFORM_TESTS
    test_slot_to_coeff
    test_slot_to_coeff_contract
    test_coeff_to_slot_contract
    test_coeff_to_slot
    test_double_hoisted_bsgs
)

set(M2424_EVALROUND_AND_BOOTSTRAP_TESTS
    test_evalround_planner
    test_canonical_embedding_reference
    test_evalround_plus_coeff_to_slot
    test_mod_raise
)

foreach(target IN LISTS
        M2424_ADAPTER_AND_SECURITY_TESTS
        M2424_LINEAR_TRANSFORM_TESTS
        M2424_EVALROUND_AND_BOOTSTRAP_TESTS)
    m2424_add_test(${target})
endforeach()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    target_link_libraries(test_coeff_to_slot PRIVATE m2424::evalmod_analysis)
    target_compile_definitions(test_coeff_to_slot PRIVATE M2424_TEST_EVALMOD_PIPELINE)
    foreach(target IN ITEMS test_evalmod_contract test_evalmod_synthesis test_evalmod_backend
                            test_evalmod_certificates test_evalmod_feasibility test_evalround_reference
                            test_evalround_backend test_evalround_certificates)
        add_executable(${target} ${target}.cpp)
        target_compile_features(${target} PRIVATE cxx_std_17)
        target_link_libraries(${target} PRIVATE m2424_test_support m2424::evalmod_analysis)
        add_test(NAME ${target} COMMAND ${target})
        set_tests_properties(${target} PROPERTIES TIMEOUT 60)
    endforeach()
    add_executable(test_evalround_k64_synthesis test_evalround_k64_synthesis.cpp)
    target_compile_features(test_evalround_k64_synthesis PRIVATE cxx_std_17)
    target_link_libraries(test_evalround_k64_synthesis PRIVATE m2424_test_support m2424::evalmod_analysis)
    add_test(NAME test_evalround_k64_synthesis COMMAND test_evalround_k64_synthesis)
    set_tests_properties(test_evalround_k64_synthesis PROPERTIES TIMEOUT 3600)
    set_tests_properties(test_evalmod_synthesis test_evalmod_backend test_evalmod_certificates
                         PROPERTIES TIMEOUT 1800)
    # These tests are deterministic but sensitive to host load. Baseline Release
    # runs exceeded the old generic 60 s timeout without failing their checks.
    set_tests_properties(test_evalmod_feasibility PROPERTIES TIMEOUT 300)
endif()

set_tests_properties(test_coeff_to_slot PROPERTIES TIMEOUT 600)

set_tests_properties(test_evalround_plus_coeff_to_slot PROPERTIES TIMEOUT 600)

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    set_tests_properties(test_evalround_backend test_evalround_certificates PROPERTIES TIMEOUT 600)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    foreach(target IN ITEMS test_slot_to_coeff test_slot_to_coeff_contract)
        target_link_libraries(${target} PRIVATE m2424::evalmod_analysis)
        target_compile_definitions(${target} PRIVATE M2424_TEST_STC_BACKEND)
        set_tests_properties(${target} PROPERTIES TIMEOUT 1800)
    endforeach()
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_coeff_to_slot_cert test_coeff_to_slot_cert.cpp)
    target_link_libraries(test_coeff_to_slot_cert PRIVATE m2424_test_support m2424::evalmod_analysis)
    add_test(NAME test_coeff_to_slot_cert COMMAND test_coeff_to_slot_cert)
    set_tests_properties(test_coeff_to_slot_cert PROPERTIES TIMEOUT 1800)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_bootstrap_end_to_end test_bootstrap_end_to_end.cpp)
    target_link_libraries(test_bootstrap_end_to_end PRIVATE m2424_test_support m2424::evalmod_analysis)
    add_test(NAME test_bootstrap_end_to_end COMMAND test_bootstrap_end_to_end)
    set_tests_properties(test_bootstrap_end_to_end PROPERTIES TIMEOUT 600)

    add_executable(bench_bootstrap_end_to_end_fixture bench_bootstrap_end_to_end_fixture.cpp)
    target_compile_features(bench_bootstrap_end_to_end_fixture PRIVATE cxx_std_17)
    target_link_libraries(bench_bootstrap_end_to_end_fixture
        PRIVATE m2424_test_support m2424::evalmod_analysis)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_sparse_bootstrap test_sparse_bootstrap.cpp)
    target_link_libraries(test_sparse_bootstrap PRIVATE m2424_test_support m2424::evalmod_analysis)
    add_test(NAME test_sparse_bootstrap COMMAND test_sparse_bootstrap)
    set_tests_properties(test_sparse_bootstrap PROPERTIES TIMEOUT 600)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_sparse_coeff_to_slot test_sparse_coeff_to_slot.cpp)
    target_link_libraries(test_sparse_coeff_to_slot PRIVATE m2424_test_support m2424::evalmod_analysis)
    add_test(NAME test_sparse_coeff_to_slot COMMAND test_sparse_coeff_to_slot)
    set_tests_properties(test_sparse_coeff_to_slot PROPERTIES TIMEOUT 2400)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_evalround_polynomial test_evalround_polynomial.cpp)
    target_link_libraries(test_evalround_polynomial PRIVATE m2424_test_support m2424::evalmod_analysis)
    add_test(NAME test_evalround_polynomial COMMAND test_evalround_polynomial)
    set_tests_properties(test_evalround_polynomial PROPERTIES TIMEOUT 600)
endif()
