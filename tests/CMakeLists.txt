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
    test_bootstrap_layout_v2
    test_bootstrap_precision_model
    test_bootstrap_scalable_refresh
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

add_executable(diagnose_bootstrap_invariants diagnose_bootstrap_invariants.cpp)
target_compile_features(diagnose_bootstrap_invariants PRIVATE cxx_std_17)
target_link_libraries(diagnose_bootstrap_invariants PRIVATE m2424)

add_executable(diagnose_dft_precision_model diagnose_dft_precision_model.cpp)
target_compile_features(diagnose_dft_precision_model PRIVATE cxx_std_17)
target_link_libraries(diagnose_dft_precision_model PRIVATE m2424)

add_executable(diagnose_rotation_noise_floor diagnose_rotation_noise_floor.cpp)
target_compile_features(diagnose_rotation_noise_floor PRIVATE cxx_std_17)
target_link_libraries(diagnose_rotation_noise_floor PRIVATE m2424)

add_executable(diagnose_bootstrap_one_cycle_precision diagnose_bootstrap_one_cycle_precision.cpp)
target_compile_features(diagnose_bootstrap_one_cycle_precision PRIVATE cxx_std_17)
target_link_libraries(diagnose_bootstrap_one_cycle_precision PRIVATE m2424)

add_executable(diagnose_bootstrap_lattice_invariant diagnose_bootstrap_lattice_invariant.cpp)
target_compile_features(diagnose_bootstrap_lattice_invariant PRIVATE cxx_std_17)
target_link_libraries(diagnose_bootstrap_lattice_invariant PRIVATE m2424)

add_executable(diagnose_stc_scaled_modup_lattice diagnose_stc_scaled_modup_lattice.cpp)
target_compile_features(diagnose_stc_scaled_modup_lattice PRIVATE cxx_std_17)
target_link_libraries(diagnose_stc_scaled_modup_lattice PRIVATE m2424)

add_executable(diagnose_bootstrap_encrypted_vs_reference diagnose_bootstrap_encrypted_vs_reference.cpp)
target_compile_features(diagnose_bootstrap_encrypted_vs_reference PRIVATE cxx_std_17)
target_link_libraries(diagnose_bootstrap_encrypted_vs_reference PRIVATE m2424)
