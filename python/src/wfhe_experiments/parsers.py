"""Parsers for benchmark output."""

import csv
import io
import re


class CsvFormatError(ValueError):
    """Raised when a benchmark did not emit a rectangular CSV table."""


def parse_csv_output(stdout: str) -> list[dict[str, str]]:
    """Parse a CSV table printed by a benchmark executable."""

    content = stdout.strip()
    if not content:
        raise CsvFormatError("benchmark produced no CSV output")

    reader = csv.DictReader(io.StringIO(content))
    if not reader.fieldnames or any(not field for field in reader.fieldnames):
        raise CsvFormatError("CSV header is missing or contains an empty field")

    records: list[dict[str, str]] = []
    for line_number, row in enumerate(reader, start=2):
        if None in row:
            raise CsvFormatError(f"CSV row {line_number} has too many columns")
        if any(value is None for value in row.values()):
            raise CsvFormatError(f"CSV row {line_number} has too few columns")
        if any(value.strip() for value in row.values()):
            records.append(dict(row))

    if not records:
        raise CsvFormatError("CSV contains a header but no data rows")
    return records


_COEFF_TO_SLOT_CANDIDATE = re.compile(
    r"^(?P<candidate>candidate_\d+),"
    r"(?P<factorization>(?:\d+-)+)"
    r"diag=(?P<estimated_diagonals>\d+),"
    r"exec_rot=(?P<estimated_rotations>\d+),"
    r"keys=(?P<estimated_keys>\d+),"
    r"actual_mul=(?P<actual_multiplications>\d+),"
    r"actual_rot=(?P<actual_rotations>\d+),"
    r"actual_keys=(?P<actual_keys>\d+)$"
)


def parse_coeff_to_slot_output(stdout: str) -> list[dict[str, str]]:
    """Normalize the multi-section output emitted by bench_coeff_to_slot."""

    records: list[dict[str, str]] = []
    section = ""
    for line_number, raw_line in enumerate(stdout.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line == "phase,ms":
            section = "phase"
            continue
        if line == "metric,value":
            section = "metric"
            continue

        candidate = _COEFF_TO_SLOT_CANDIDATE.fullmatch(line)
        if candidate:
            row = candidate.groupdict()
            row["factorization"] = row["factorization"].rstrip("-")
            records.append({"record_type": "candidate", **row})
            continue

        if section in {"phase", "metric"}:
            parts = line.split(",", maxsplit=1)
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise CsvFormatError(f"invalid CoeffToSlot row {line_number}: {line}")
            measurement, value = parts
            if section == "phase":
                unit = "ms"
            elif measurement in {"first_error", "second_error"}:
                unit = "absolute_error"
            elif measurement == "prepared_bytes":
                unit = "bytes"
            else:
                unit = "count"
            records.append(
                {
                    "record_type": section,
                    "measurement": measurement,
                    "value": value,
                    "unit": unit,
                }
            )
            continue

        raise CsvFormatError(f"unknown CoeffToSlot row {line_number}: {line}")

    if not records:
        raise CsvFormatError("benchmark produced no CoeffToSlot records")
    if not any(row.get("record_type") == "phase" for row in records):
        raise CsvFormatError("CoeffToSlot output has no timing section")
    return records
