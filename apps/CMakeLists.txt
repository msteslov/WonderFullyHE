option(M2424_BUILD_APPS "Build m2424 command-line demos and measurement tools" ON)
option(M2424_BUILD_RESEARCH_APPS "Build experimental/research bootstrap harnesses" OFF)

if(NOT M2424_BUILD_APPS)
    return()
endif()

function(m2424_add_app target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endfunction()

set(M2424_STABLE_APPS
    bench_chain_accuracy
    bench_ckks
    bench_parameter_planner
    bench_parallel_throughput
    demo_abft
    demo_basic
    demo_checked_pipeline
    demo_client_compute_roundtrip
    demo_eval_mod_polynomial
    demo_galois_key_optimization
    demo_mod_raise
    demo_noise_growth
    demo_precision_profiles
    demo_profile_report
    demo_secure_stats
    demo_security_report
)

set(M2424_RESEARCH_APPS
    bench_bootstrap_full
    bench_bootstrap_layout_physical_gate
    bench_bootstrap_layout_planner
    bench_bootstrap_multi_cycle
    bench_bootstrap_one_case
    bench_bootstrap_parts
    bench_bootstrap_period_model
    bench_bootstrap_planner
    bench_bootstrap_prescaled_coeff_to_slot
    bench_bootstrap_profile_budget
    bench_bootstrap_reference_path
    bench_bootstrap_refresh
    bench_bootstrap_scale_strategy
    bench_bootstrap_scaling
    bench_bootstrap_slots_to_coeffs_first
    bench_bootstrap_stc_budget
    bench_bootstrap_stc_deep_probe
    bench_bootstrap_validation
    demo_bootstrap_cipher_path
    demo_bootstrap_diagonals
    demo_bootstrap_end_to_end
    demo_bootstrap_pipeline
    demo_bootstrap_prototype
)

foreach(target IN LISTS M2424_STABLE_APPS)
    m2424_add_app(${target})
endforeach()

if(M2424_BUILD_RESEARCH_APPS)
    foreach(target IN LISTS M2424_RESEARCH_APPS)
        m2424_add_app(${target})
    endforeach()
endif()
