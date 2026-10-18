import pytest

from wfhe_experiments.parsers import (
    CsvFormatError,
    parse_coeff_to_slot_output,
    parse_csv_output,
)


def test_parse_csv_output():
    records = parse_csv_output("operation,time_ms\nencode,1.25\ndecode,2.5\n")

    assert records == [
        {"operation": "encode", "time_ms": "1.25"},
        {"operation": "decode", "time_ms": "2.5"},
    ]


@pytest.mark.parametrize("content", ["", "a,b\n", "a,b\n1,2,3\n"])
def test_invalid_csv_is_rejected(content):
    with pytest.raises(CsvFormatError):
        parse_csv_output(content)


def test_parse_coeff_to_slot_output():
    records = parse_coeff_to_slot_output(
        """phase,ms
context,48.5
apply_warm,300.0
metric,value
prepared_bytes,1024
first_error,8e-11
candidate_0,6-5-5-3-diag=266,exec_rot=81,keys=71,actual_mul=476,actual_rot=134,actual_keys=55
"""
    )

    assert records[0] == {
        "record_type": "phase",
        "measurement": "context",
        "value": "48.5",
        "unit": "ms",
    }
    assert records[2]["unit"] == "bytes"
    assert records[-1]["factorization"] == "6-5-5-3"
    assert records[-1]["actual_rotations"] == "134"
