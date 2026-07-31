function(m2424_add_test target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424)
    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES TIMEOUT 60)
endfunction()

set(M2424_TESTS
    test_accuracy
    test_adapter_failures
    test_bootstrap_candidates
    test_canonical_embedding_reference
    test_coeff_to_slot_contract
    test_coeff_to_slot
    test_ckks_primitives
    test_mod_raise
)

foreach(target IN LISTS M2424_TESTS)
    m2424_add_test(${target})
endforeach()

set_tests_properties(test_coeff_to_slot PROPERTIES TIMEOUT 600)
