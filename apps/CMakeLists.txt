option(M2424_BUILD_APPS "Build m2424 command-line demos and measurement tools" ON)

if(NOT M2424_BUILD_APPS)
    return()
endif()

function(m2424_add_app target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424)
endfunction()

set(M2424_BENCHMARKS
    bench_chain_accuracy
    bench_bootstrap_candidates
    bench_ckks
    bench_coeff_to_slot
    bench_parallel_throughput
)

set(M2424_DEMOS
    demo_abft
    demo_basic
    demo_noise_growth
    demo_precision_profiles
    demo_profile_report
    demo_secure_stats
    demo_security_report
)

foreach(target IN LISTS M2424_BENCHMARKS M2424_DEMOS)
    m2424_add_app(${target})
endforeach()

set(M2424_REPORTS
    report_evalmod_synthesis
    report_evalmod_feasibility
    report_evalround_k64_synthesis
)

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    foreach(target IN LISTS M2424_REPORTS)
        add_executable(${target} ${target}.cpp)
        target_compile_features(${target} PRIVATE cxx_std_17)
        target_link_libraries(${target} PRIVATE m2424::evalmod_analysis)
    endforeach()
endif()
