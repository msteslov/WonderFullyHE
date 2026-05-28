add_executable(test_smoke test_smoke.cpp)
target_compile_features(test_smoke PRIVATE cxx_std_17)
target_link_libraries(test_smoke PRIVATE m2424)

add_executable(test_accuracy test_accuracy.cpp)
target_compile_features(test_accuracy PRIVATE cxx_std_17)
target_link_libraries(test_accuracy PRIVATE m2424)

add_executable(test_adapter_failures test_adapter_failures.cpp)
target_compile_features(test_adapter_failures PRIVATE cxx_std_17)
target_link_libraries(test_adapter_failures PRIVATE m2424)

add_test(NAME test_smoke COMMAND test_smoke)
add_test(NAME test_accuracy COMMAND test_accuracy)
add_test(NAME test_adapter_failures COMMAND test_adapter_failures)
