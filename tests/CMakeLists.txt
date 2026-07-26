function(m2424_add_test target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424)
    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES TIMEOUT 60)
endfunction()

function(m2424_add_research_test target)
    m2424_add_test(${target})
    set_tests_properties(${target} PROPERTIES LABELS research TIMEOUT 300)
endfunction()

set(M2424_TESTS
    test_accuracy
    test_adapter_failures
    test_bootstrap_dft
    test_bootstrap_experiment
    test_bootstrap_layout_v2
    test_bootstrap_precision_model
    test_bootstrap_scale_design
    test_bootstrap_stc_reference
    test_mod1_approximation
    test_mod1_circuit
    test_operation_budget
    test_parameter_planner
    test_smoke
)

foreach(target IN LISTS M2424_TESTS)
    m2424_add_test(${target})
endforeach()

option(M2424_BUILD_RESEARCH_TESTS "Build experimental bootstrap tests and diagnostics" OFF)

if(M2424_BUILD_RESEARCH_TESTS)
    m2424_add_research_test(test_bootstrap_scalable_refresh)
    m2424_add_research_test(test_bootstrap_multi_cycle_precision)

    set(M2424_RESEARCH_DIAGNOSTICS
        diagnose_bootstrap_encrypted_vs_reference
        diagnose_bootstrap_invariants
        diagnose_bootstrap_lattice_invariant
        diagnose_bootstrap_one_cycle_precision
        diagnose_dft_precision_model
        diagnose_rotation_noise_floor
        diagnose_stc_scaled_modup_lattice
        diagnose_wide_mod1_plain
    )

    foreach(target IN LISTS M2424_RESEARCH_DIAGNOSTICS)
        add_executable(${target} ${target}.cpp)
        target_compile_features(${target} PRIVATE cxx_std_17)
        target_link_libraries(${target} PRIVATE m2424)
    endforeach()
endif()
