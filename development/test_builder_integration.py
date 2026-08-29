"""
Exercises the real jeolData.append_specdata_via_builder() and
createJsonDict_via_builder() methods added 2026-08-27 — not
reimplementations of their logic, the actual bound methods on the actual
class, via jeolData.__new__() to skip __init__'s HDF5/RDKit file
processing (which needs a real .jjh5 file we don't have here) and set
just the attributes those two methods actually read.

Includes a direct regression check for the DEPT135/HSQC_CLIPCOSY case
that used to crash (KeyError) or silently vanish before the 2026-08-27
dialog + counter-dict fix — this test would have failed against the
original code.
"""

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, str(Path(__file__).parent))

import pytest
import simpleNMRjeolTools_v8 as v8


def _scalar_env(value, datatype="string"):
    return {"datatype": datatype, "count": 1, "data": {"0": value}}


def _string_list_env(strings, datatype="special"):
    return {"datatype": datatype, "count": len(strings), "data": {str(i): s for i, s in enumerate(strings)}}


def _indexed_env(records, datatype="dataframe"):
    return {"datatype": datatype, "count": len(records), "data": {str(i): r for i, r in enumerate(records)}}


def _minimal_dataset(expt_type: str) -> dict:
    """A synthetic dataset_with_peaks entry with just enough structure for
    append_specdata_via_builder to build a valid spectrum block."""
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


def _make_stand_in_jeol_data(expt_types: list[str]) -> "v8.jeolData":
    """
    Build a jeolData instance without running __init__ (which needs a real
    .jjh5 file), populating only the attributes append_specdata_via_builder
    and createJsonDict_via_builder actually read.
    """
    obj = v8.jeolData.__new__(v8.jeolData)

    datasets_with_peaks = {}
    chosen_entries = []
    for i, expt_type in enumerate(expt_types):
        key = f"expt_{i}"
        datasets_with_peaks[key] = _minimal_dataset(expt_type)
        # matches the real "<...metadata...> <key> <type>" shape closely
        # enough to exercise the same split-based parsing append_specdata
        # uses (expt_id_parts[-2] = key, expt_id_parts[-1] = type)
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


@pytest.mark.skipif(not v8._SIMPLENMR_BUILDER_AVAILABLE, reason="simplenmr_builder not installed")
def test_dept135_and_hsqc_clipcosy_no_longer_crash_or_vanish():
    """
    Direct regression test for the confirmed 2026-08-27 bug: before the
    fix, DEPT135 raised KeyError in the old expts_type_count dict lookup,
    and HSQC_CLIPCOSY (as the old mislabeled "HSQCCLIPCOSY") never
    matched a real NMREXPERIMENTS token and silently vanished. Both must
    now appear as real blocks in the built payload.
    """
    obj = _make_stand_in_jeol_data(["HSQC", "DEPT135", "HSQC_CLIPCOSY"])
    builder = v8.SimpleNMRBuilder(source="jeol")

    obj.append_specdata_via_builder(builder)  # must not raise

    blocks = builder.spectra.build_blocks()
    assert "DEPT135_0" in blocks
    assert "HSQC_CLIPCOSY_0" in blocks
    assert "HSQC_0" in blocks


@pytest.mark.skipif(not v8._SIMPLENMR_BUILDER_AVAILABLE, reason="simplenmr_builder not installed")
def test_unrecognized_token_raises_clear_error_not_keyerror():
    """
    A token that ISN'T a real NMREXPERIMENTS type should fail loudly and
    specifically at add_block() time (ValueError), never as a bare
    KeyError deep inside dict indexing, and never silently.
    """
    obj = _make_stand_in_jeol_data(["HSQC", "SOME_MADE_UP_TYPE"])
    builder = v8.SimpleNMRBuilder(source="jeol")

    with pytest.raises(ValueError, match="not a recognized NMREXPERIMENTS token"):
        obj.append_specdata_via_builder(builder)


@pytest.mark.skipif(not v8._SIMPLENMR_BUILDER_AVAILABLE, reason="simplenmr_builder not installed")
def test_create_json_dict_via_builder_full_roundtrip():
    """
    End-to-end: the real createJsonDict_via_builder() method, called on a
    stand-in jeolData instance, produces a payload that (a) omits
    nmrAssignments/exptIdentifiers as required for jeol, (b) contains the
    HSQC block, and (c) passes the real vendored validator.
    """
    obj = _make_stand_in_jeol_data(["HSQC", "COSY"])

    payload = obj.createJsonDict_via_builder()

    assert "nmrAssignments" not in payload
    assert "exptIdentifiers" not in payload
    assert "HSQC_0" in payload
    assert "COSY_0" in payload


@pytest.mark.skipif(not v8._SIMPLENMR_BUILDER_AVAILABLE, reason="simplenmr_builder not installed")
def test_working_directory_and_filename_survive_the_builder():
    """
    Regression test for a confirmed 2026-08-27 bug: createJsonDict_via_builder()
    originally omitted workingDirectory/workingFilename entirely. The
    server itself doesn't need them (server_consumed: false), but
    submit_to_server()'s local post-submission code
    (_write_html_result()/_log_submission_error()) reads them straight
    back OUT of the payload after a successful server round trip to know
    where to write the HTML result and error log — omitting them crashed
    that local step with a bare KeyError even on a successful submission,
    surfacing as "Error submitting to simpleNMR server: 'workingFilename'"
    despite the server having already responded with status 200.
    """
    obj = _make_stand_in_jeol_data(["HSQC"])

    payload = obj.createJsonDict_via_builder()

    assert payload["workingDirectory"]["data"]["0"] == "/fake/working/dir"
    assert payload["workingFilename"]["data"]["0"] == "fake_sample"
