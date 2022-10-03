from CPExpectedDistance.p_expected_distance import expected_distance
import RNA
import pytest


TESTSEQS = [
    "AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUUACGAAAGGCUGUAAAAUCAAUUAUUCACCACAGGGGGCCCCCGUGUCUAG",
    "AAUACGUACGCAUUGCACAUCG"
]


@pytest.mark.parametrize(
    "seq,temperature,constrains,expected",
    [
        (TESTSEQS[0], 37, None, 10.182507219601332),
        (TESTSEQS[0], 15, None, 10.099126824664843),
        (TESTSEQS[0], 37, [50, 65], 14.686194628267502),

    ]
)
def test_expected_distance(seq, temperature, constrains, expected):
    md = RNA.md(temperature=temperature)
    fc = RNA.fold_compound(seq, md)
    if constrains is not None:
        for x in range(constrains[0], constrains[1]):
            fc.hc_add_up(x)
    arr = expected_distance(fc)
    assert arr[0, -1] == expected
