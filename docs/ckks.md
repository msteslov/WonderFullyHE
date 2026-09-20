# Базовые операции CKKS

`SealAdapter` предоставляет интерфейс CKKS без непосредственного использования
типов Microsoft SEAL. Для подключения используется заголовок
`include/m2424/m2424.hpp` и цель `m2424::m2424`.

## Минимальный пример

Фрагмент вычисляет покомпонентный квадрат зашифрованного вектора:

```cpp
#include "m2424/m2424.hpp"

auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
adapter.generateKeys(true, true);

auto encrypted = adapter.encrypt(adapter.encode({1.0, 2.0, 3.0}));
auto product = adapter.multiply(encrypted, encrypted);
auto reduced = adapter.relinearize(product);
auto squared = adapter.rescaleToNext(reduced);
auto decoded = adapter.decode(adapter.decrypt(squared));
```

Первые три значения результата приближённо равны `1`, `4`, `9`.
Расшифровка здесь используется для получения результата; её не требуется
выполнять на стороне гомоморфного вычисления.
Полные примеры находятся в demo_basic.cpp (`apps/demo_basic.cpp`)
и demo_secure_stats.cpp (`apps/demo_secure_stats.cpp`).

## Операции и совместимость

| Группа | Методы |
|---|---|
| Кодирование | `encode`, `encodeFor`, `encodeScalarFor`, `decode` |
| Шифрование | `encrypt`, `decrypt` |
| Аддитивные операции | `add`, `sub`, `addPlain`, `subPlain` |
| Умножение | `multiply`, `multiplyPlain`, `relinearize`, `rescaleToNext` |
| Перестановки и уровни | `rotate`, `modSwitchTo` |

Имена API используют `camelCase`. Суффикс `For` означает кодирование
plaintext на уровне и с масштабом заданного шифротекста.
Пользователь обязан учитывать совместимость контекста, уровня и масштаба,
наличие ключей и запас modulus chain. Операции CKKS не обеспечивают
неограниченную глубину вычислений.

Сериализация позволяет разделить контексты шифрования, вычисления с evaluation
keys и расшифровки. Набор ключей должен соответствовать исполняемым операциям.

## Профили и точность

`m2424::profiles` содержит `fast_demo_ckks`, `basic_ckks`, `balanced_ckks`,
`depth_ckks`, `high_precision_ckks`. Фактические параметры доступны через:

```bash
./build/bin/demo_profile_report
./build/bin/demo_security_report
```

Сравнение с эталоном выполняется средствами `accuracy`:
`max_abs_error`, `mean_abs_error`, `compare(expected, actual, tolerance)`.
Точность определяется профилем, диапазоном входа и последовательностью операций.
Тестовый допуск не является универсальной гарантией для произвольной программы.

Модуль `abft` реализует checksum-проверки сложения, вычитания, покомпонентного
умножения и ротации. Они проверяют численную согласованность и не заменяют
криптографическую аутентификацию результата.

Команды сборки и перечень измерительных приложений приведены в
[инструкции по сборке](build.md).
