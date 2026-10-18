# Python experiments

Этот пакет запускает C++-бенчмарки WonderFullyHE, сохраняет сырые данные и
метаданные воспроизводимости, строит статистическую сводку и графики.

## Установка и проверка

Из корня репозитория:

```bash
python3 -m venv python/.venv
python/.venv/bin/pip install -e 'python[dev]'
python/.venv/bin/ruff check python
python/.venv/bin/pytest python
```

Перед запуском экспериментов необходимо собрать C++-цели. Быстрая проверка:

```bash
python/.venv/bin/wfhe-experiments run python/configs/smoke.json
```

Полный baseline с тремя CKKS-профилями и кандидатами bootstrap:

```bash
python/.venv/bin/wfhe-experiments run python/configs/baseline.json
```

Исследовательский прогон с точностью по глубине, parallel throughput и пятью
факторизациями CoeffToSlot:

```bash
cmake -S . -B build-python-research -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DM2424_ENABLE_EVALMOD_ANALYSIS=ON
cmake --build build-python-research -j
python/.venv/bin/wfhe-experiments run python/configs/research.json
```

End-to-end benchmark использует явно обозначенный test-only fixture. Он проверяет
и измеряет полный арифметический маршрут, но не является свидетельством безопасности
или production-сертификатом lift.

Исследовательский вопрос, гипотезы, метрики и ограничения зафиксированы в
[`docs/research_protocol.md`](../docs/research_protocol.md).

Каждый запуск создаёт отдельный каталог в `python/results/` с файлами
`manifest.json`, `records.csv`, `records.json`, сырыми stdout/stderr, отчётом и
графиками. Manifest фиксирует commit, наличие незакоммиченных изменений, команды,
версии Python и платформы.
