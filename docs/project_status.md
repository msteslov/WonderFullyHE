# Статус проекта

WonderFullyHE — учебно-исследовательский прототип библиотеки для защищённых вычислений на CKKS поверх Microsoft SEAL.

## Реализовано

- Адаптер `SealAdapter`, скрывающий низкоуровневые типы Microsoft SEAL.
- CKKS-профили `basic_ckks` и `depth_ckks`.
- Операции `encode`, `encrypt`, `decrypt`, `decode`.
- Гомоморфные операции `add`, `sub`, `mul_relin_rescale`, `rotate`.
- Метрики точности `max_abs_error`, `mean_abs_error`, `tolerance`.
- ABFT-проверки для `add`, `sub`, `mul` и `rotate`.
- Демонстрация защищённой статистической агрегации.
- Анализ расхода вычислительной глубины до остановки на `x^32`.
- Bootstrapping-конвейер `ModRaise -> CoeffToSlot -> EvalMod -> SlotToCoeff` с отчётом по `scale`, `chain_index` и размеру ciphertext.
- Benchmark времени операций, ошибок и размеров ciphertext/ключей.
- Security report по CKKS-профилям относительно лимитов Microsoft SEAL `tc128`, `tc192` и `tc256`.
- CMake-сборка, smoke-тесты, GitLab CI и GitHub Actions CI.

## В работе

- Вычислительные блоки `CoeffToSlot`, `EvalMod`, `SlotToCoeff`.
- End-to-end refresh ciphertext с проверками `Dec(c') ~= Dec(c)` и `level(c') > level(c)`.
- Разделение клиентского и серверного контекстов для облачного сценария.
- Безопасная сериализация ciphertext/ключей и проверка совместимости параметров.
- Расширение ABFT на цепочки операций.

## Ограничения

- Проект не является production-ready криптографической библиотекой.
- ABFT-модуль контролирует вычислительную согласованность, но не заменяет криптографическую аутентификацию результата.
- Для реального облачного сценария требуется отдельное управление ключами, формат обмена данными и анализ side-channel рисков.
