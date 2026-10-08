"""
Ethylbenzene, HSQC + HMBC + COSY only (real export, 2026-10-04).

The quaternary carbon C2 has no HSQC peak, so simpleNMR reports it as
the *mean* 13C coordinate of its three HMBC peaks (144.394-144.401 ->
144.398). Exact matching found nothing; the long-range step must pick
those three peaks and only those.
"""

import json
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

import beautifuljason as bjason
from simplenmrjeol import assignment_writer as aw

DATA = Path(__file__).parent / "data" / "ethylbenzene_py.json"


@pytest.fixture
def snmr():
    return json.loads(DATA.read_text())


def _document(snmr):
    old = snmr["oldjsondata"]

    def spectrum(key, name):
        nuc = tuple(old[key]["nucleus"])
        peaks = [NS(pos=(p["delta2"], p["delta1"]), h5_group=NS(name=f"{name}/{i}"))
                 for i, p in enumerate(old[key]["peaks"]["data"].values())]
        return NS(ndim=2, spec_info=NS(nuclides=nuc), peaks=peaks)

    graph = snmr["molgraph"]
    atoms = [NS(type=bjason.Molecule.Atom.Type(n["atom_number"]), bonded_atom_numbers=[]) for n in graph["nodes"]]
    for l in graph["links"]:
        atoms[l["source"]].bonded_atom_numbers.append(l["target"])
        atoms[l["target"]].bonded_atom_numbers.append(l["source"])

    shift = lambda a, v: NS(nums=[a], mark=None, ignored_auto=False, ignored_user=False,
                            get_value_error_pair=lambda _m, v=v: (v, 1.0))
    # one shift per atom, as JASON stored it (14 assignments before the fix = 7 C + 7 H)
    c = [shift(a, v) for a, v in [(0, 127.7), (1, 144.0), (2, 127.7), (3, 128.3), (4, 125.6),
                                  (5, 128.3), (6, 29.4), (7, 15.7)]]
    h = [shift(a, v) for a, v in [(0, 7.2), (2, 7.2), (3, 7.3), (4, 7.2), (5, 7.3), (6, 2.6), (7, 1.2)]]
    molecule = NS(atoms=atoms, spectra=[NS(nucleus=aw.C13, shifts=c), NS(nucleus=aw.H1, shifts=h)])
    return NS(mol_data=[molecule],
              nmr_data=[spectrum("COSY_0", "cosy"), spectrum("HMBC_0", "hmbc"), spectrum("HSQC_0", "hsqc")])


def test_quaternary_carbon_from_hmbc_mean(snmr):
    _, tuples, result = aw.build_assignments(_document(snmr), snmr)
    assert result.problems == []
    assert result.n_assignments == len(tuples) == 15          # 14 before + C2
    # one shift -> one signal: C2 on the HMBC peak closest to simpleNMR's value
    c2 = [t for t in tuples if t[0].nums == [1]]
    assert len(c2) == 1
    assert c2[0][1].h5_group.name.startswith("hmbc/") and c2[0][2] == 1
    assert round(c2[0][1].pos[1], 3) == 144.4
    row = [r for r in result.rows if r.atom == "C2"][0]
    assert "mean of 3 = simpleNMR value" in row.status


def test_symmetric_partners(snmr):
    nodes, _, _ = aw.parse_simplenmr_export(snmr)
    labels = [n["label"] for n in nodes]
    assert "C1=C3" in labels and "C6=C4" in labels


def test_unused_1d_13c_in_document_does_not_block_hmbc(snmr):
    """The JASON document may hold a 1D 13C spectrum simpleNMR didn't use."""
    doc = _document(snmr)
    unrelated = [NS(pos=(x,), h5_group=NS(name=f"c13/{i}")) for i, x in enumerate([77.16, 128.0, 143.9])]
    doc.nmr_data.append(NS(ndim=1, spec_info=NS(nuclides=("13C",)), peaks=unrelated))
    _, tuples, result = aw.build_assignments(doc, snmr)
    assert result.problems == [] and len(tuples) == 15

