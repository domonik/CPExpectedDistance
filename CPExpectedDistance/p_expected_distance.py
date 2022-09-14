from CPExpectedDistance.c_expected_distance import cp_expected_distance
from typing import Dict


def expected_distance(seq: str, md_config: Dict = None):
    if md_config is None:
        md_config = {}
    return cp_expected_distance(seq, md_config)
