option(M2424_BUILD_APPS "Build m2424 command-line demos and measurement tools" ON)

if(NOT M2424_BUILD_APPS)
    return()
endif()

function(m2424_add_app target)
    add_executable(${target} ${target}.cpp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_link_libraries(${target} PRIVATE m2424)
endfunction()

set(M2424_STABLE_APPS
    bench_chain_accuracy
    bench_bootstrap_candidates
    bench_ckks
    bench_coeff_to_slot
    bench_parallel_throughput
    demo_abft
    demo_basic
    demo_noise_growth
    demo_precision_profiles
    demo_profile_report
    demo_secure_stats
    demo_security_report
)

foreach(target IN LISTS M2424_STABLE_APPS)
    m2424_add_app(${target})
endforeach()

if(M2424_ENABLE_EVALMOD_ANALYSIS)
    add_executable(report_evalmod_synthesis report_evalmod_synthesis.cpp)
    target_compile_features(report_evalmod_synthesis PRIVATE cxx_std_17)
    target_link_libraries(report_evalmod_synthesis PRIVATE m2424::evalmod_analysis)
    add_executable(report_evalmod_feasibility report_evalmod_feasibility.cpp)
    target_compile_features(report_evalmod_feasibility PRIVATE cxx_std_17)
    target_link_libraries(report_evalmod_feasibility PRIVATE m2424::evalmod_analysis)
endif()
