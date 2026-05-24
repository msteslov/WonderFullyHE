# Статус проекта

WonderFullyHE — учебно-исследовательский прототип библиотеки для защищённых вычислений на CKKS поверх Microsoft SEAL.

## Реализовано

- Адаптер `SealAdapter`, скрывающий низкоуровневые типы Microsoft SEAL.
- Готовые CKKS-профили в `m2424::profiles`: `fast_demo_ckks`, `basic_ckks`, `balanced_ckks`, `depth_ckks`, `high_precision_ckks`.
- Операции `encode`, `encrypt`, `decrypt`, `decode`.
- Гомоморфные операции `add`, `sub`, `mul_relin_rescale`, `rotate`.
- Plaintext-операции `add_plain`, `sub_plain`, `mul_plain_rescale`.
- Кодирование plaintext на уровне ciphertext через `encode_like` и `encode_scalar_like`.
- Генерация только нужных Galois-ключей для заданных rotation steps.
- Выравнивание ciphertext по уровню и масштабу через `mod_switch_to` и `match_level_and_scale`.
- Метрики точности `max_abs_error`, `mean_abs_error`, `tolerance`.
- ABFT-проверки для `add`, `sub`, `mul` и `rotate`.
- Rotation-based `LinearTransform`.
- `sum_slots` для суммы слотов в первом CKKS-слоте.
- `PolynomialEvaluator` для вычисления полиномов от ciphertext.
- Демонстрация защищённой статистической агрегации.
- Демонстрация оптимизации Galois-ключей: полный набор против ограниченного набора ротаций.
- Анализ расхода вычислительной глубины до остановки на `x^32`.
- Bootstrapping-конвейер `ModRaise -> CoeffToSlot -> EvalMod -> SlotToCoeff` с отчётом по `scale`, `chain_index` и размеру ciphertext.
- Benchmark времени операций, ошибок и размеров ciphertext/ключей.
- Benchmark строительных блоков bootstrapping: plaintext multiplication, linear transform, slot summation, polynomial evaluation.
- Benchmark параллельной обработки независимых ciphertext с разделением setup/runtime.
- Security report по CKKS-профилям относительно лимитов Microsoft SEAL `tc128`, `tc192` и `tc256`.
- CMake alias target `m2424::m2424` для подключения библиотеки через `add_subdirectory`.
- CMake-сборка, smoke-тесты, GitLab CI и GitHub Actions CI.

## В работе

- Подстановка рассчитанных матриц в `CoeffToSlot` и `SlotToCoeff`.
- Подстановка коэффициентов полинома в `EvalMod`.
- End-to-end refresh ciphertext с проверками `Dec(c') ~= Dec(c)` и `level(c') > level(c)`.
- Разделение клиентского и серверного контекстов для облачного сценария.
- Безопасная сериализация ciphertext/ключей и проверка совместимости параметров.
- Расширение ABFT на цепочки операций.

## Ограничения

- ABFT-модуль контролирует вычислительную согласованность, но не заменяет криптографическую аутентификацию результата.
- Для реального облачного сценария требуется отдельное управление ключами, формат обмена данными и анализ side-channel рисков.
