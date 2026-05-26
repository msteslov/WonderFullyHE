# Статус проекта

WonderFullyHE — учебно-исследовательский прототип библиотеки для защищённых вычислений на CKKS поверх Microsoft SEAL.

## Реализовано

- Адаптер `SealAdapter`, скрывающий низкоуровневые типы Microsoft SEAL.
- Готовые CKKS-профили в `m2424::profiles`: `fast_demo_ckks`, `basic_ckks`, `balanced_ckks`, `depth_ckks`, `high_precision_ckks`.
- Операции `encode`, `encrypt`, `decrypt`, `decode`.
- Гомоморфные операции `add`, `sub`, `mul_relin_rescale`, `rotate`.
- Plaintext-операции `add_plain`, `sub_plain`, `mul_plain`, `mul_plain_rescale`.
- Кодирование plaintext на уровне ciphertext через `encode_like`, `encode_scalar_like` и `encode_scalar_at_scale_like`.
- Генерация только нужных Galois-ключей для заданных rotation steps.
- Сериализация и загрузка публичного ключа, секретного ключа, Relin/Galois-ключей и ciphertext.
- Разделённый roundtrip-сценарий: шифрование данных, вычисление без secret key, расшифрование результата.
- Выравнивание ciphertext по уровню и масштабу через `mod_switch_to` и `match_level_and_scale`.
- Метрики точности `max_abs_error`, `mean_abs_error`, `tolerance`.
- ABFT-проверки для `add`, `sub`, `mul` и `rotate`.
- Rotation-based `LinearTransform`.
- `sum_slots` для суммы слотов в первом CKKS-слоте.
- `PolynomialEvaluator` для вычисления полиномов от ciphertext.
- Демонстрация защищённой статистической агрегации.
- Демонстрация оптимизации Galois-ключей: полный набор против ограниченного набора ротаций.
- Анализ расхода вычислительной глубины до остановки на `x^32`.
- Низкоуровневый CKKS `ModRaise` для расширения ciphertext к первой RNS-базе.
- Bootstrapping-конвейер `ModRaise -> CoeffToSlot -> eval_mod_normalization -> EvalMod -> SlotToCoeff -> post_refresh_mod_raise` с отчётом по `scale`, `chain_index` и размеру ciphertext.
- Публичный API `Bootstrapper::refresh(cipher, slots, tolerance)` и `Bootstrapper::refresh_rotation_steps(slots)`.
- Нормализация входа `EvalMod` по амплитуде после `CoeffToSlot`.
- End-to-end demo, где после refresh выполняется следующая plaintext-операция с rescale.
- Benchmark публичного refresh-пути.
- Benchmark времени операций, ошибок и размеров ciphertext/ключей.
- Benchmark строительных блоков bootstrapping: plaintext multiplication, linear transform, slot summation, polynomial evaluation.
- Benchmark параллельной обработки независимых ciphertext с разделением setup/runtime.
- Security report по CKKS-профилям относительно лимитов Microsoft SEAL `tc128`, `tc192` и `tc256`.
- CMake alias target `m2424::m2424` для подключения библиотеки через `add_subdirectory`.
- CMake-сборка, smoke-тесты, GitLab CI и GitHub Actions CI.

## В работе

- Расширение ABFT на цепочки операций.
- Sweep-benchmark для refresh на разных размерах входа и параметрах.

## Ограничения

- ABFT-модуль контролирует вычислительную согласованность, но не заменяет криптографическую аутентификацию результата.
- Для реального внешнего применения требуется отдельное управление ключами, формат обмена данными и анализ side-channel рисков.
