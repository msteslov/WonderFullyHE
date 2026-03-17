add_executable(test_smoke test_smoke.cpp)
target_compile_features(test_smoke PRIVATE cxx_std_17)
target_link_libraries(test_smoke PRIVATE m2424)

add_test(NAME test_smoke COMMAND test_smoke)
