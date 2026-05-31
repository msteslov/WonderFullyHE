function(m2424_add_test target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424)
    add_test(NAME ${target} COMMAND ${target})
endfunction()

set(M2424_TESTS
    test_accuracy
    test_adapter_failures
    test_bootstrap_dft
    test_bootstrap_scalable_refresh
    test_bootstrap_scale_design
    test_mod1_circuit
    test_operation_budget
    test_parameter_planner
    test_smoke
)

foreach(target IN LISTS M2424_TESTS)
    m2424_add_test(${target})
endforeach()
