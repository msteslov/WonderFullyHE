function(m2424_add_test target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424)
    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES TIMEOUT 60)
endfunction()

set(M2424_TESTS
    test_sparse_security
    test_slot_to_coeff
    test_slot_to_coeff_contract
    test_accuracy
    test_adapter_failures
    test_bootstrap_candidates
    test_bootstrap_contract
    test_evalround_planner
    test_canonical_embedding_reference
    test_coeff_to_slot_contract
    test_coeff_to_slot
    test_evalround_plus_coeff_to_slot
    test_ckks_primitives
    test_double_hoisted_bsgs
    test_mod_raise
)

foreach(target IN LISTS M2424_TESTS)
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
        target_link_libraries(${target} PRIVATE m2424::evalmod_analysis)
        add_test(NAME ${target} COMMAND ${target})
        set_tests_properties(${target} PROPERTIES TIMEOUT 60)
    endforeach()
    set_tests_properties(test_evalmod_synthesis test_evalmod_backend test_evalmod_certificates
                         PROPERTIES TIMEOUT 1800)
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
    target_link_libraries(test_coeff_to_slot_cert PRIVATE m2424::evalmod_analysis)
    add_test(NAME test_coeff_to_slot_cert COMMAND test_coeff_to_slot_cert)
    set_tests_properties(test_coeff_to_slot_cert PROPERTIES TIMEOUT 1800)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_bootstrap_end_to_end test_bootstrap_end_to_end.cpp)
    target_link_libraries(test_bootstrap_end_to_end PRIVATE m2424::evalmod_analysis)
    add_test(NAME test_bootstrap_end_to_end COMMAND test_bootstrap_end_to_end)
    set_tests_properties(test_bootstrap_end_to_end PROPERTIES TIMEOUT 600)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_sparse_bootstrap test_sparse_bootstrap.cpp)
    target_link_libraries(test_sparse_bootstrap PRIVATE m2424::evalmod_analysis)
    add_test(NAME test_sparse_bootstrap COMMAND test_sparse_bootstrap)
    set_tests_properties(test_sparse_bootstrap PROPERTIES TIMEOUT 600)
endif()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(test_sparse_coeff_to_slot test_sparse_coeff_to_slot.cpp)
    target_link_libraries(test_sparse_coeff_to_slot PRIVATE m2424::evalmod_analysis)
    add_test(NAME test_sparse_coeff_to_slot COMMAND test_sparse_coeff_to_slot)
    set_tests_properties(test_sparse_coeff_to_slot PROPERTIES TIMEOUT 600)
endif()
