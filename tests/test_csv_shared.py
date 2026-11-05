"""Shared CSV text decoding and read-only array conversion."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import experimental_csv, rydberg_experiment, rydberg_experimental_csv  # noqa: E402
from gabes import vapor_cell  # noqa: E402
from gabes._arrays import readonly_float_array  # noqa: E402
from gabes._csv_numeric import decode_csv_text  # noqa: E402

TEXT = "Detuning [MHz],Signal\n-1,0.2\n0,0.8\n1,0.2\n"


@pytest.mark.parametrize("data, encoding", [
    (TEXT.encode("utf-16"), "utf-16"),
    (b"\xfe\xff" + TEXT.encode("utf-16-be"), "utf-16"),
    (TEXT.encode("utf-8-sig"), "utf-8-sig"),
    (TEXT.encode("utf-8"), "utf-8-sig"),
    (("메모\n" + TEXT).encode("cp949"), "cp949"),
    (b"\xff\x81 note\n" + TEXT.encode("ascii"), "latin-1"),
])
def test_decode_csv_text_branches(data, encoding):
    text, used = decode_csv_text(data)
    assert used == encoding
    assert text.endswith(TEXT) and not text.startswith("﻿")


def test_each_loader_keeps_its_own_error_for_malformed_utf16():
    bad = b"\xff\xfe" + TEXT.encode("utf-16-le")[:-1]
    with pytest.raises(KeyError, match="Malformed UTF-16"):
        decode_csv_text(bad, KeyError)
    with pytest.raises(experimental_csv.ExperimentalCSVError, match="Malformed UTF-16"):
        experimental_csv.load_experimental_csv(bad)
    with pytest.raises(rydberg_experimental_csv.RydbergCSVError, match="Malformed UTF-16"):
        rydberg_experimental_csv.load_eit_csv(bad)


def test_readonly_float_array_copy_flag():
    source = np.arange(4.0)
    copied = readonly_float_array(source)
    assert not copied.flags.writeable and source.flags.writeable
    assert not np.shares_memory(copied, source)
    frozen = readonly_float_array(source, copy=False)
    assert frozen is source and not source.flags.writeable
    converted = readonly_float_array([1, 2, 3], copy=False)
    assert converted.dtype == np.float64 and not converted.flags.writeable


def test_callers_use_the_shared_converter():
    assert rydberg_experimental_csv._readonly is readonly_float_array
    assert vapor_cell._readonly is readonly_float_array
    assert rydberg_experiment._readonly is readonly_float_array
    assert experimental_csv._readonly.func is readonly_float_array
    assert experimental_csv._readonly.keywords == {"copy": False}
