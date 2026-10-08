"""
assignment_writer on a molecule with symmetry and no 1D 13C spectrum
(4-fluorobenzotrifluoride, real run 2026-10-01). This case failed in
the first version: every 13C value came from the HSQC/HMBC 13C axis and
was looked for only in 1D 13C spectra, and the symmetric carbons C6/C7
were exported only through their partners C11/C10.
"""

import json
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

import beautifuljason as bjason
from simplenmrjeol import assignment_writer as aw

DATA = Path(__file__).parent / "data" / "fluoro_cf3_py.json"


@pytest.fixture
def snmr():
    return json.loads(DATA.read_text())


def _spectrum(nuclides, positions, name):
    peaks = [NS(pos=pos, h5_group=NS(name=f"{name}/{i}")) for i, pos in enumerate(positions)]
    return NS(ndim=len(nuclides), spec_info=NS(nuclides=nuclides), peaks=peaks)


def _shift(nums, value):
    return NS(nums=nums, mark=None, ignored_auto=False, ignored_user=False,
              get_value_error_pair=lambda _m, v=value: (v, 1.0))


def _document(snmr, combine_symmetric):
    old = snmr["oldjsondata"]
    pk = lambda key: old[key]["peaks"]["data"].values()
    hsqc = _spectrum(("1H", "13C"), [(p["delta2"], p["delta1"]) for p in pk("HSQC_0")], "hsqc")
    hmbc = _spectrum(("1H", "13C"), [(p["delta2"], p["delta1"]) for p in pk("HMBC_0")], "hmbc")
    cosy = _spectrum(("1H", "1H"), [(p["delta2"], p["delta1"]) for p in pk("COSY_0")], "cosy")

    graph = snmr["molgraph"]
    atoms = [NS(type=bjason.Molecule.Atom.Type(n["atom_number"]), bonded_atom_numbers=[]) for n in graph["nodes"]]
    for l in graph["links"]:
        atoms[l["source"]].bonded_atom_numbers.append(l["target"])
        atoms[l["target"]].bonded_atom_numbers.append(l["source"])

    if combine_symmetric:     # JASON lists equivalent atoms in one shift
        c = [_shift([1], 124.4), _shift([4], 125.9), _shift([5, 10], 127.5), _shift([6, 9], 116.2), _shift([7], 165.0)]
        h = [_shift([5, 10], 7.6), _shift([6, 9], 7.2)]
    else:                     # one shift per atom
        c = [_shift([a], v) for a, v in [(1, 124.4), (4, 125.9), (5, 127.5), (6, 116.2), (7, 165.0), (9, 116.2), (10, 127.5)]]
        h = [_shift([a], v) for a, v in [(5, 7.6), (6, 7.2), (9, 7.2), (10, 7.6)]]
    molecule = NS(atoms=atoms, spectra=[NS(nucleus=aw.C13, shifts=c), NS(nucleus=aw.H1, shifts=h)])
    return NS(mol_data=[molecule], nmr_data=[cosy, hsqc, hmbc])


def test_symmetric_partners_are_folded_in(snmr):
    nodes, _, _ = aw.parse_simplenmr_export(snmr)
    assert [(n["label"], n["atom_ids"]) for n in nodes] == [
        ("C2", [1]), ("C5", [4]), ("C8", [7]), ("C10=C7", [9, 6]), ("C11=C6", [10, 5])]


@pytest.mark.parametrize("combine, expected", [(False, 11), (True, 7)])
def test_all_carbons_assigned_without_1d_13c(snmr, combine, expected):
    _, tuples, result = aw.build_assignments(_document(snmr, combine), snmr)
    assert result.problems == []
    assert result.n_assignments == len(tuples) == expected
    assert result.n_carbons == 7 and result.n_symmetric == 2
    # protonated carbons use the HSQC 13C axis (dim 1); quaternary carbons the HMBC 13C axis
    for t in tuples:
        name = t[1].h5_group.name
        assert not name.startswith("cosy/")
        if name.startswith("hmbc/"):
            assert t[2] == 1
    quaternary = {r.atom: r.status for r in result.rows if r.nucleus == "13C" and "long-range" in r.status}
    assert set(quaternary) == {"C2", "C5", "C8"}
    assert "mean of 2" in quaternary["C8"]          # two HMBC peaks share 164.685 ppm; one is used


def test_long_range_can_be_switched_off(snmr, monkeypatch):
    monkeypatch.setattr(aw, "ASSIGN_LONG_RANGE_CARBONS", False)
    _, _, result = aw.build_assignments(_document(snmr, False), snmr)
    assert sorted(r.atom for r in result.problems) == ["C2", "C5", "C8"]
