"""
Tests for assignment_writer.build_assignments() against a mock JASON
document built from a real simpleNMR export (alpha-ionone). The real
.jjh5 (~146 MB) is too large to commit; the mock reproduces what
build_assignments() reads: molecule atoms/bonds, predicted shifts, and
the 13C 1D / HSQC / HMBC peak tables. The peak positions come from the
export's own oldjsondata, which matched the real document exactly
(max difference 0) when checked on 2026-09-30.
"""

import json
import os
import time
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

import beautifuljason as bjason
from simplenmrjeol.assignment_writer import (
    AtomMismatchError,
    C13,
    H1,
    build_assignments,
    find_export,
)

DATA = Path(__file__).parent / "data" / "alpha_ionone_py.json"


@pytest.fixture
def snmr():
    return json.loads(DATA.read_text())


def _spectrum(nuclides, positions, name):
    peaks = [NS(pos=pos, h5_group=NS(name=f"{name}/{i}")) for i, pos in enumerate(positions)]
    return NS(ndim=len(nuclides), spec_info=NS(nuclides=nuclides), peaks=peaks)


def _shift(atom, value, mark=None):
    return NS(nums=[atom], mark=mark, ignored_auto=False, ignored_user=False,
              get_value_error_pair=lambda _method, v=value: (v, 1.0))


def _document(snmr):
    old = snmr["oldjsondata"]
    peaks = lambda key: old[key]["peaks"]["data"].values()
    c13 = _spectrum(("13C",), [(p["delta1"],) for p in peaks("C13_1D_0")], "c13")
    hsqc = _spectrum(("1H", "13C"), [(p["delta2"], p["delta1"]) for p in peaks("HSQC_0")], "hsqc")
    hmbc = _spectrum(("1H", "13C"), [(p["delta2"], p["delta1"]) for p in peaks("HMBC_0")], "hmbc")

    graph = snmr["molgraph"]
    atoms = [NS(type=bjason.Molecule.Atom.Type(n["atom_number"]), bonded_atom_numbers=[])
             for n in graph["nodes"]]
    for link in graph["links"]:
        atoms[link["source"]].bonded_atom_numbers.append(link["target"])
        atoms[link["target"]].bonded_atom_numbers.append(link["source"])

    c_shifts, h_shifts = [], []
    for node in snmr["nodes_now"]:
        if node["symbol"] != "C":
            continue
        i = node["id"]
        c_shifts.append(_shift(i, node["ppm_calculated"]))
        if node["numProtons"] == 2:          # JASON marks CH2 protons dn/up
            h_shifts += [_shift(i, 1.30, "dn"), _shift(i, 1.35, "up")]
        elif node["numProtons"]:
            h_shifts.append(_shift(i, 1.0))

    molecule = NS(atoms=atoms, spectra=[NS(nucleus=C13, shifts=c_shifts), NS(nucleus=H1, shifts=h_shifts)])
    return NS(mol_data=[molecule], nmr_data=[hmbc, c13, hsqc])


def test_all_values_assigned(snmr):
    doc = _document(snmr)
    _, tuples, result = build_assignments(doc, snmr)
    assert result.n_assignments == len(tuples) == 25
    assert result.problems == []
    # 13C tuples point at the 13C 1D peaks; 1H tuples at HSQC peaks, never HMBC
    for t in tuples:
        name = t[1].h5_group.name
        if len(t) == 2:
            assert name.startswith("c13/")
        else:
            assert name.startswith("hsqc/") and t[2] == 0


def test_ch2_pairing(snmr):
    _, _, result = build_assignments(_document(snmr), snmr)
    c8 = [(r.shift, r.simplenmr) for r in result.rows if r.atom == "C8" and r.nucleus == "1H"]
    assert c8 == [("H8(dn)", pytest.approx(1.2076, abs=1e-4)), ("H8(up)", pytest.approx(1.4388, abs=1e-4))]
    c12 = [r.shift for r in result.rows if r.atom == "C12" and r.nucleus == "1H"]
    assert c12 == ["H12(dn)", "H12(up)"]       # one value, used for both shifts


def test_moved_peak_is_reported_not_bridged(snmr):
    doc = _document(snmr)
    c13 = doc.nmr_data[1]
    target = [n["ppm"] for n in snmr["nodes_now"] if n["id"] == 11][0]
    for p in c13.peaks:
        if abs(p.pos[0] - target) < 1e-6:
            p.pos = (p.pos[0] + 0.02,)
    _, tuples, result = build_assignments(doc, snmr)
    assert len(tuples) == 24
    assert [(r.atom, r.status) for r in result.problems] == [("C12", "no exact 13C peak")]


def test_atom_mismatch_stops(snmr):
    doc = _document(snmr)
    doc.mol_data[0].atoms[13].type = bjason.Molecule.Atom.Type.N
    with pytest.raises(AtomMismatchError):
        build_assignments(doc, snmr)


def test_find_export_ignores_stale_files(tmp_path):
    old = tmp_path / "input_py.json"
    old.write_text("{}")
    past = time.time() - 60
    os.utime(old, (past, past))
    since = time.time() - 1
    assert find_export(tmp_path, since) is None
    new = tmp_path / "input_py.json"
    new.write_text("{}")
    assert find_export(tmp_path, since) == new
    assert find_export(tmp_path / "missing", since) is None



def one_signal_per_shift(tuples):
    signals = {}
    for t in tuples:
        signals.setdefault(id(t[0]), set()).add((id(t[1]), t[2] if len(t) > 2 else None))
    return all(len(v) == 1 for v in signals.values())


def test_each_shift_has_one_signal(snmr):
    """JASON: 'A chemical shift cannot be assigned to more than one signal'."""
    _, tuples, _ = build_assignments(_document(snmr), snmr)
    assert one_signal_per_shift(tuples)
