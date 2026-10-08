"""
Route 1 of assignment_writer: assignments and peak moves from the export's
peak_map (simpleNMR server with JEOL IDs).

Fixtures are a real ethylbenzene run (October 2026):
  * ethylbenzene_jeol_input.json  - the input sent to simpleNMR: every peak's
    JEOL ID and original position (used to build a mock JASON document);
  * ethylbenzene_peakmap_py.json  - the viewer's export, with peak_map.
The expected values are the ATOMS and MOVES tables that gave the correct
display in JASON (simpleNMR_JASON_findings_ethylbenzene.ipynb).
"""

import copy
import json
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
import pytest

import beautifuljason as bjason
from simplenmrjeol import assignment_writer as aw

DATA = Path(__file__).parent / "data"
EXPORT = json.loads((DATA / "ethylbenzene_peakmap_py.json").read_text())
JEOL_INPUT = json.loads((DATA / "ethylbenzene_jeol_input.json").read_text())

# atom (1-based) -> (1D 13C peak ID, HSQC peak ID), confirmed in JASON
ATOMS = {
    1: ("{981e8a05-cc8d-479a-a660-57d1de40ed90}", "{b02fad61-a90c-43ae-9115-ebbc664d9306}"),  # ortho
    2: ("{db5d5e8d-a444-4c9f-8652-052226a52d1a}", None),                                      # quaternary
    3: ("{981e8a05-cc8d-479a-a660-57d1de40ed90}", "{b02fad61-a90c-43ae-9115-ebbc664d9306}"),  # ortho
    4: ("{e0ce5f38-c916-452d-ac68-55a58caf9cc6}", "{fc3bd16a-e67c-47d6-bc1d-73daf2d26bc1}"),  # meta
    5: ("{7dd23078-a81f-4aca-a9e8-ceef2c501e24}", "{f80dacb4-d950-47eb-8ad9-7a8b46596427}"),  # para
    6: ("{e0ce5f38-c916-452d-ac68-55a58caf9cc6}", "{fc3bd16a-e67c-47d6-bc1d-73daf2d26bc1}"),  # meta
    7: ("{60d8daa1-0725-4220-9bc5-0ed42dc2e803}", "{bfb0db5f-1680-4585-aa52-955107815816}"),  # CH2
    8: ("{84fce28e-910a-4c13-bf36-b2839f6ad45a}", "{49dad087-f74e-45da-ada3-1bd45abb1bbc}"),  # CH3
}

# JEOL ID -> (Pos[0], Pos[1]) after moving, confirmed in JASON
MOVES = {
    "{1774e688-a307-4f0c-92ad-b60b6da1c87b}": (7.4012, 144.3054),   # HMBC meta-H -> C2
    "{bb84d936-9adf-44bd-8c46-ebd0a7aad069}": (1.3676, 144.3054),   # HMBC CH3 -> C2
    "{17aa3cdf-e8f8-478b-b07c-124e637c561c}": (1.3676, 28.9897),    # HMBC CH3 -> C7
    "{9cb4a1c8-fb58-4864-a418-a446a0803e6e}": (2.7775, 144.3054),   # HMBC CH2 -> C2
    "{82292e56-7132-4af1-a2df-269c3383b84b}": (2.7775, 127.9292),   # HMBC CH2 -> ortho
    "{00a416e8-712d-48b2-b645-73536f816579}": (7.2956, 127.9292),   # HMBC para-H -> ortho
    "{c722cbc6-e190-4a44-b790-757bcb2691ab}": (2.7775, 1.3676),     # COSY CH2 - CH3
    "{0ce9e917-3c5f-4d11-8297-93b930ab33fd}": (7.4012, 7.3226),     # COSY meta - ortho
    "{e48c637d-1420-4c76-941c-39c6c30dd5a8}": (7.4012, 7.2956),     # COSY meta - para
}


def _document():
    """Mock JASON document: real peak IDs/positions, one predicted shift per atom."""
    spectra = []
    for block in JEOL_INPUT.values():
        if not (isinstance(block, dict) and block.get("datatype") == "nmrspectrum"):
            continue
        nuclides = tuple([block["nucleus"]] if isinstance(block["nucleus"], str) else block["nucleus"])
        peaks = []
        for p in block["peaks"]["data"].values():
            # 1D: the ppm is in delta1 (json_converter swaps 1D deltas); 2D: Pos = (delta2, delta1)
            pos = np.array([p["delta1"], 0.0, 0.0] if len(nuclides) == 1 else [p["delta2"], p["delta1"], 0.0])
            peaks.append(NS(id=p["jason_id"], pos=pos, h5_group=NS(attrs={"Pos": pos})))
        spectra.append(NS(spec_info=NS(nuclides=nuclides), peaks=peaks))

    graph = EXPORT["molgraph"]
    atoms = [NS(type=bjason.Molecule.Atom.Type(n["atom_number"]), bonded_atom_numbers=[]) for n in graph["nodes"]]
    for l in graph["links"]:
        atoms[l["source"]].bonded_atom_numbers.append(l["target"])
        atoms[l["target"]].bonded_atom_numbers.append(l["source"])

    shift = lambda atom, v: NS(nums=[atom - 1], mark=None, ignored_auto=False, ignored_user=False,
                               get_value_error_pair=lambda _m, v=v: (v, 1.0))
    c = [shift(a, v) for a, v in [(1, 127.7), (2, 144.0), (3, 127.7), (4, 128.3), (5, 125.6),
                                  (6, 128.3), (7, 29.4), (8, 15.7)]]
    h = [shift(a, v) for a, v in [(1, 7.1), (3, 7.1), (4, 7.1), (5, 7.1), (6, 7.1), (7, 2.6), (8, 1.2)]]
    molecule = NS(atoms=atoms, spectra=[NS(nucleus=aw.C13, shifts=c), NS(nucleus=aw.H1, shifts=h)])
    return NS(mol_data=[molecule], nmr_data=spectra)


def _assigned(tuples):
    """{(nucleus, atom 1-based): (peak ID, dim)} from the assignment tuples."""
    out = {}
    for t in tuples:
        nucleus = "C" if any(t[0] is s for s in _assigned.c_shifts) else "H"
        out[(nucleus, t[0].nums[0] + 1)] = (t[1].id, t[2] if len(t) > 2 else None)
    return out


def _build(export):
    doc = _document()
    _assigned.c_shifts = doc.mol_data[0].spectra[0].shifts
    return aw.build_assignments(doc, export)


def test_assignments_by_id():
    _, tuples, result = _build(EXPORT)
    assert result.problems == []
    got = _assigned(tuples)
    for atom, (c13_id, hsqc_id) in ATOMS.items():
        assert got[("C", atom)] == (c13_id, None)
        if hsqc_id:
            assert got[("H", atom)] == (hsqc_id, 0)       # HSQC 1H axis
    assert len(tuples) == 15


def test_moves_match_the_recipe():
    _, _, result = _build(EXPORT)
    got = {m.jason_id: (round(m.new[0], 4), round(m.new[1], 4)) for m in result.moves}
    assert got == MOVES
    assert "9 HMBC/COSY peak(s) will be moved" in result.summary_text()


def _swap_in_viewer(export, atom_a, atom_b):
    """What the viewer's commitMoves does: data (ppm, H1_ppm, id) changes atoms; atomNumber,
    symmetry and coordinates stay with the atom."""
    export = copy.deepcopy(export)
    nodes = {n["atomNumber"]: n for n in export["nodes_now"]}
    a, b = nodes[atom_a], nodes[atom_b]
    for key in ("ppm", "H1_ppm", "id"):
        a[key], b[key] = b[key], a[key]
    return export


def test_shift_swapped_in_viewer_follows_atom_number():
    """Ortho (atom 1, symmetric with 3) swapped with para (atom 5): the para peaks go to
    atoms 1 and 3, the ortho peaks to atom 5. (The old writer used the node id as the
    atom and put them back, with para data on one ortho atom as well.)"""
    _, tuples, result = _build(_swap_in_viewer(EXPORT, 1, 5))
    assert result.problems == []
    got = _assigned(tuples)
    para_c, para_h = ATOMS[5]
    ortho_c, ortho_h = ATOMS[1]
    assert got[("C", 1)] == got[("C", 3)] == (para_c, None)
    assert got[("H", 1)] == got[("H", 3)] == (para_h, 0)
    assert got[("C", 5)] == (ortho_c, None)
    assert got[("H", 5)] == (ortho_h, 0)
    # correlation peaks move to the same positions: they belong to the signals, not the atoms
    assert {m.jason_id: (round(m.new[0], 4), round(m.new[1], 4)) for m in result.moves} == MOVES


def test_missing_id_is_reported_not_guessed():
    export = copy.deepcopy(EXPORT)
    export["peak_map"]["correlations"][0]["jason_id"] = "{00000000-0000-0000-0000-000000000000}"
    _, _, result = _build(export)
    assert len(result.moves) == 8
    assert [r.status for r in result.problems] == ["JEOL ID {00000000-0000-0000-0000-000000000000} not found in the document"]


def test_write_moves_peaks_and_logs(tmp_path, monkeypatch):
    doc = _document()
    _assigned.c_shifts = doc.mol_data[0].spectra[0].shifts
    written = []
    doc.mol_data[0].set_assignments = written.append

    class FakeDocument:
        def __init__(self, path, mode="r"):
            pass

        def __enter__(self):
            return doc

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(aw.bjason, "Document", FakeDocument)
    monkeypatch.setattr(aw.bjason.Molecule, "set_assignments", lambda *a: None, raising=False)
    jjh5 = tmp_path / "input.jjh5"
    jjh5.write_bytes(b"")

    result = aw.write_assignments(jjh5, EXPORT)

    assert len(written[0]) == 15
    by_id = {p.id: p for s in doc.nmr_data for p in s.peaks}
    for jason_id, (pos0, pos1) in MOVES.items():
        pos = by_id[jason_id].h5_group.attrs["Pos"]
        assert (round(pos[0], 4), round(pos[1], 4)) == (pos0, pos1)
    log = json.loads((tmp_path / "input_peak_moves.json").read_text())
    assert len(log) == 9 and {e["jason_id"] for e in log} == set(MOVES)
    meta_c2 = next(e for e in log if e["jason_id"] == "{1774e688-a307-4f0c-92ad-b60b6da1c87b}")
    assert [round(x, 4) for x in meta_c2["original"][:2]] == [7.4039, 144.4006]
