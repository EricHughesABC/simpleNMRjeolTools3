"""
Exercises the real jeolData.append_specdata() and createJsonDict()
methods — not reimplementations of their logic, the actual bound methods
on the actual class, via jeolData.__new__() to skip __init__'s HDF5/RDKit
file processing (which needs a real .jjh5 file we don't have here) and
set just the attributes those two methods actually read.

Includes a direct regression check for the DEPT135/HSQC_CLIPCOSY case
that used to crash (KeyError) or silently vanish before the 2026-08-27
dialog + counter-dict fix — this test would have failed against the
original code.

simplenmr_builder is a hard dependency of this package (see
json_converter.py) — there's no skip marker for "not installed" the way
earlier versions of this test had, since that's no longer a supported
state; if it's missing, collection itself fails with a clear ImportError
pointing at the fix, same as running the real program would.
"""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from simplenmr_builder import SimpleNMRBuilder
from simplenmrjeol.json_converter import jeolData


def _scalar_env(value, datatype="string"):
    return {"datatype": datatype, "count": 1, "data": {"0": value}}


def _string_list_env(strings, datatype="special"):
    return {"datatype": datatype, "count": len(strings), "data": {str(i): s for i, s in enumerate(strings)}}


def _indexed_env(records, datatype="dataframe"):
    return {"datatype": datatype, "count": len(records), "data": {str(i): r for i, r in enumerate(records)}}


def _minimal_dataset(expt_type: str) -> dict:
    """A synthetic dataset_with_peaks entry with just enough structure for
    append_specdata to build a valid spectrum block."""
    return {
        "spec_info": {
            "pulsesequence": "fake_pulseseq",
            "experimenttype": expt_type,
            "datafilename": "fake.jjh5",
            "solvent": "CDCl3",
            "expt_fn": "fake_expt",
            "specfrequency": [400.0],
            "type": "2D" if expt_type not in ("H1_1D", "C13_1D") else "1D",
            "temperature": 298.0,
            "nucleus": ["1H", "13C"],
        },
        "peaks": {"datatype": "peaks", "count": 0, "data": {}},
        "integrals": {"datatype": "integrals", "count": 0, "data": {}, "normValue": 1},
        "multiplets": {"datatype": "multiplets", "count": 0, "data": {}, "normValue": 1},
    }


def _make_stand_in_jeol_data(expt_types):
    """
    Build a jeolData instance without running __init__ (which needs a real
    .jjh5 file), populating only the attributes append_specdata() and
    createJsonDict() actually read.
    """
    obj = jeolData.__new__(jeolData)

    datasets_with_peaks = {}
    chosen_entries = []
    for i, expt_type in enumerate(expt_types):
        key = f"expt_{i}"
        datasets_with_peaks[key] = _minimal_dataset(expt_type)
        chosen_entries.append(f"[1H, 13C] fakepulseseq {key} {expt_type}")

    obj.datasets_with_peaks = datasets_with_peaks
    obj.chosenSpectra = _string_list_env(chosen_entries, datatype="chosenSpectra")
    obj.spectraWithPeaks = _string_list_env(chosen_entries, datatype="spectraWithPeaks")

    obj.smiles = _scalar_env("CCO")
    obj.molfile = _scalar_env("fake molfile contents")
    obj.hostname = _scalar_env("TESTHOST-1234")
    obj.MNOVAcalcMethod = _scalar_env("JEOL Predict")
    obj.carbonCalcPositionsMethod = _scalar_env("Calculated Positions")
    obj.workingDirectory = _scalar_env("/fake/working/dir", datatype="workingDirectory")
    obj.workingFilename = _scalar_env("fake_sample", datatype="workingFilename")
    obj.simulated_annealing = True
    obj.ml_consent = False

    obj.allatomsInfo = _indexed_env([{"atom_idx": 0, "atomNumber": 1, "numProtons": 3}])
    obj.carbonatomsInfo = _indexed_env([{"atom_idx": 0, "atomNumber": 1}])
    obj.c13predictions = _indexed_env([])

    return obj


def test_dept135_and_hsqc_clipcosy_no_longer_crash_or_vanish():
    obj = _make_stand_in_jeol_data(["HSQC", "DEPT135", "HSQC_CLIPCOSY"])
    builder = SimpleNMRBuilder(source="jeol")

    obj.append_specdata(builder)  # must not raise

    blocks = builder.spectra.build_blocks()
    assert "DEPT135_0" in blocks
    assert "HSQC_CLIPCOSY_0" in blocks
    assert "HSQC_0" in blocks


def test_unrecognized_token_raises_clear_error_not_keyerror():
    obj = _make_stand_in_jeol_data(["HSQC", "SOME_MADE_UP_TYPE"])
    builder = SimpleNMRBuilder(source="jeol")

    with pytest.raises(ValueError, match="not a recognized NMREXPERIMENTS token"):
        obj.append_specdata(builder)


def test_create_json_dict_full_roundtrip():
    obj = _make_stand_in_jeol_data(["HSQC", "COSY"])

    payload = obj.createJsonDict()

    assert "nmrAssignments" not in payload
    assert "exptIdentifiers" not in payload
    assert "HSQC_0" in payload
    assert "COSY_0" in payload


def test_working_directory_and_filename_survive_the_builder():
    obj = _make_stand_in_jeol_data(["HSQC"])

    payload = obj.createJsonDict()

    assert payload["workingDirectory"]["data"]["0"] == "/fake/working/dir"
    assert payload["workingFilename"]["data"]["0"] == "fake_sample"
