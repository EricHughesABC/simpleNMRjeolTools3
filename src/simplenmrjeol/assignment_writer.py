"""
assignment_writer.py

Writes the assignments exported from the simpleNMR viewer (the
``<workingFilename>_py.json`` file) back into a JASON ``.jjh5`` document
using ``Molecule.set_assignments()`` (beautifuljason >= 1.3.0b1).

Two routes:

1. **By JEOL ID** (exports from a simpleNMR server that writes ``peak_map``).
   The server records, for every carbon group, the JEOL IDs of its 1D 13C
   peak and HSQC peak(s), and for every HMBC/COSY peak it used, the
   position simpleNMR snapped it to. The writer then
     * assigns each carbon to its 1D 13C peak (else its HSQC peak's 13C
       axis, else one of its HMBC peaks' 13C axis) and each proton to its
       HSQC peak's 1H axis - all looked up by ID, no ppm matching;
     * moves each HMBC/COSY peak to its snapped position, so JASON groups
       it with the assigned signals and draws the correlation (the recipe
       confirmed by hand on ethylbenzene, October 2026). The original
       positions are logged next to the .jjh5.
   The atom comes from the node's ``atomNumber``; its ``id`` names the
   carbon group and travels with it when the user moves a shift in the
   viewer.

2. **By ppm matching** (older exports without ``peak_map``), described
   below. Everything in the export was originally read from the same
   ``.jjh5``, so nothing is approximated:

* Atoms: the node's ``atomNumber`` - 1 is the JASON atom. Elements and
  heavy-atom bonds are compared and an ``AtomMismatchError`` is raised if
  they differ; there is no guessing.
* 13C: each node's ``ppm`` must equal a peak position, looked for in
  this order:
    1. a 1D 13C spectrum -> ``(shift, peak)``;
    2. the 13C axis of a 1H/13C 2D peak whose 1H coordinate is also one
       of the carbon's protons (the HSQC peak) -> ``(shift, peak, c_dim)``;
    3. the 13C axis of the long-range (HMBC-type) 1H/13C peaks that
       belong to this carbon -> ``(shift, peak, c_dim)``, if
       ``ASSIGN_LONG_RANGE_CARBONS`` and simpleNMR did not use a 1D 13C
       spectrum (read from the export's ``oldjsondata``; the document may
       still contain one). This is how quaternary carbons are found from HSQC/HMBC
       data alone.
  Steps 1 and 2 are exact. Step 3 cannot be: simpleNMR reports the
  *mean* 13C coordinate of a quaternary carbon's HMBC peaks (ethylbenzene
  C2: three peaks at 144.394-144.401, reported as 144.398). The peaks
  used are those within ``LONG_RANGE_TOL`` of the value and closer to it
  than to any other carbon's value; the report says whether their mean
  reproduces the simpleNMR value exactly.
* Symmetry: simpleNMR exports one node for each set of equivalent
  carbons (``sym_atom_idx`` names the partner). The node's values are
  applied to the partner atom too, unless the partner is exported with
  its own values (e.g. diastereotopic gem-dimethyls).
* 1H: each value in ``H1_ppm`` must equal the 1H coordinate of a peak in
  a 1H/13C 2D spectrum. If several peaks match, the one whose 13C
  coordinate is closest to the carbon is used (the HSQC peak)
  -> ``(shift, peak, h_dim)``.
* CH2: predicted shifts marked dn/up are paired with the two proton
  values in ascending ppm; a single value is used for both shifts.

"Equal" means within ``PPM_EPS``, a floating-point allowance, not a
chemical tolerance. A value with no exact match is reported, never
bridged.

JASON's assignment model: each predicted shift is assigned to exactly
one signal (``set_assignments`` raises "A chemical shift cannot be
assigned to more than one signal" otherwise). Several shifts may share
one signal (symmetric atoms, CH2 protons with one value). JASON draws
HMBC/COSY correlations itself, from its merged 1D/2D signals, so those
peaks are not assigned here. (Assigning them explicitly was tried on
2026-10-05 and rejected by JASON for exactly this reason.)
"""

from __future__ import annotations

import json
import math
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import beautifuljason as bjason

PPM_EPS = 1e-4

# Assign quaternary-carbon shifts to the 13C axis of HMBC-type peaks when
# no 1D 13C peak matches. Set False if JASON rejects long-range assignments.
ASSIGN_LONG_RANGE_CARBONS = True
LONG_RANGE_TOL = 0.5      # ppm window for HMBC-type peaks around a carbon found only via long-range peaks


H1 = bjason.Molecule.Atom.NuclType.H1
C13 = bjason.Molecule.Atom.NuclType.C13
ELEMENT_H = int(bjason.Molecule.Atom.Type.H)


class AtomMismatchError(ValueError):
    """The simpleNMR molecule and the JASON molecule differ atom-for-atom."""


@dataclass
class ReportRow:
    atom: str                    # e.g. "C8"
    nucleus: str                 # "13C" or "1H"
    shift: str                   # JASON shift label, e.g. "H8(dn)"
    predicted: Optional[float]   # JASON predicted value
    simplenmr: Optional[float]   # experimental value from simpleNMR
    source: str                  # e.g. "nmr_data[4]"
    status: str                  # "assigned" or a problem description

    @property
    def ok(self) -> bool:
        return self.status.startswith("assigned")


@dataclass
class PeakMove:
    """An HMBC/COSY peak to move so JASON draws its correlation."""
    kind: str                    # "HMBC" or "COSY"
    jason_id: str
    peak: object                 # the document's peak (valid while the file is open)
    original: list               # Pos before
    new: list                    # Pos after

    def log_entry(self) -> dict:
        return {"kind": self.kind, "jason_id": self.jason_id,
                "original": [float(x) for x in self.original],
                "moved_to": [float(x) for x in self.new]}


@dataclass
class AssignmentResult:
    n_assignments: int
    n_carbons: int
    rows: list = field(default_factory=list)
    n_symmetric: int = 0          # carbons covered through a symmetric partner
    moves: list = field(default_factory=list)   # PeakMove, ID route only

    @property
    def problems(self) -> list:
        return [r for r in self.rows if not r.ok]

    def summary_text(self, max_problems: int = 15) -> str:
        carbons = f"{self.n_carbons} carbons"
        if self.n_symmetric:
            carbons += f" ({self.n_symmetric} by symmetry)"
        lines = [f"{self.n_assignments} assignments for {carbons}."]
        if self.moves:
            lines.append(f"{len(self.moves)} HMBC/COSY peak(s) will be moved to simpleNMR's positions.")
        long_range = sorted({r.atom for r in self.rows if r.ok and "long-range" in r.status})
        if long_range:
            lines.append(f"13C from HMBC-type peaks (no 1D 13C match): {', '.join(long_range)}")
            for r in self.rows:
                if r.ok and "long-range" in r.status and "simpleNMR " in r.status:
                    lines.append(f"  check {r.atom}: {r.status[len('assigned '):]}")
        problems = self.problems
        if not problems:
            lines.append("Every simpleNMR value matched a peak in the JASON document.")
        else:
            lines.append(f"{len(problems)} value(s) could not be assigned:")
            for r in problems[:max_problems]:
                value = "" if r.simplenmr is None else f" {r.simplenmr:.4f}"
                lines.append(f"  {r.atom} {r.nucleus}{value}: {r.status}")
            if len(problems) > max_problems:
                lines.append(f"  ... and {len(problems) - max_problems} more")
        return "\n".join(lines)


# ---------------------------------------------------------------- simpleNMR export

def _ppm_or_none(value, limit=1000.0) -> Optional[float]:
    """float, or None for missing values and simpleNMR's 10000 placeholder."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) and abs(value) < limit else None


def _sym_ids(value) -> list:
    """sym_atom_idx may be an int, a numeric string, '' or a comma/space separated list."""
    if value is None or value == "":
        return []
    if isinstance(value, (int, float)):
        return [int(value)]
    return [int(v) for v in str(value).replace(",", " ").split() if v.strip().lstrip("-").isdigit()]


def parse_simplenmr_export(snmr: dict):
    """Return (carbon nodes, {id: atomic number}, {frozenset bond}) from an export dict.

    Each node carries ``id`` (its atom, from ``atomNumber``), ``node_id``
    (the carbon group's identity, which travels with the data when the
    user moves a shift in the viewer, and keys ``peak_map``) and
    ``atom_ids``: its atom plus any symmetric partner that is not exported
    with values of its own.
    """
    raw = []
    for n in snmr["nodes_now"]:
        if n.get("symbol") != "C":
            continue
        h_ppms = sorted(v for v in (_ppm_or_none(x, 100) for x in n.get("H1_ppm") or []) if v is not None)
        raw.append({
            "id": int(n["atomNumber"]) - 1,   # the atom it sits on now
            "node_id": str(int(n["id"])),     # the carbon group (moves with the data)
            "c_ppm": _ppm_or_none(n.get("ppm")),
            "h_ppms": h_ppms,
            "sym": _sym_ids(n.get("sym_atom_idx")),
        })
    has_values = {n["id"] for n in raw if n["c_ppm"] is not None or n["h_ppms"]}
    covered = {s for n in raw if n["id"] in has_values for s in n["sym"] if s not in has_values}

    nodes = []
    for n in raw:
        if n["id"] in covered:          # placeholder for a symmetric partner: handled by that node
            continue
        partners = [s for s in n["sym"] if s in covered]
        label = f"C{n['id'] + 1}" + "".join(f"=C{s + 1}" for s in partners)
        nodes.append({"id": n["id"], "node_id": n["node_id"], "atom_ids": [n["id"], *partners],
                      "label": label, "c_ppm": n["c_ppm"], "h_ppms": n["h_ppms"]})
    nodes.sort(key=lambda n: n["id"])
    elements = {int(n["id"]): int(n["atom_number"]) for n in snmr["molgraph"]["nodes"]}
    bonds = {frozenset((int(l["source"]), int(l["target"]))) for l in snmr["molgraph"]["links"]}
    return nodes, elements, bonds


def find_export(html_dir: Path, since: float) -> Optional[Path]:
    """Newest ``*_py.json`` in html_dir written at or after ``since`` (a time.time() value).

    The viewer writes a fresh file on every Export press; the time check
    stops a file left over from an earlier run being used.
    """
    html_dir = Path(html_dir)
    if not html_dir.is_dir():
        return None
    fresh = [p for p in html_dir.glob("*_py.json") if p.stat().st_mtime >= since]
    return max(fresh, key=lambda p: p.stat().st_mtime) if fresh else None


def load_export(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


# ---------------------------------------------------------------- JASON document

def _check_atoms(molecule, elements: dict, bonds: set) -> None:
    j_elements = {i: int(a.type) for i, a in enumerate(molecule.atoms) if int(a.type) != ELEMENT_H}
    j_bonds = {frozenset((i, int(j))) for i in j_elements
               for j in molecule.atoms[i].bonded_atom_numbers if int(j) in j_elements}
    if j_elements != elements or j_bonds != bonds:
        raise AtomMismatchError(
            "The simpleNMR molecule and the JASON molecule do not match atom-for-atom.\n"
            f"  element differences (index, atomic number): "
            f"{sorted(set(j_elements.items()) ^ set(elements.items()))}\n"
            f"  bonds only in JASON: {sorted(tuple(sorted(b)) for b in j_bonds - bonds)}\n"
            f"  bonds only in simpleNMR: {sorted(tuple(sorted(b)) for b in bonds - j_bonds)}"
        )


def _all_peaks(document, ndim: int, nuclei: tuple) -> list:
    """[(nmr_data index, peak)] for every spectrum with this dimensionality and nuclei."""
    out = []
    for i, spectrum in enumerate(document.nmr_data):
        if spectrum.ndim == ndim and set(spectrum.spec_info.nuclides) == set(nuclei):
            try:
                out += [(i, p) for p in spectrum.peaks]
            except KeyError:        # spectrum has no peak table
                pass
    return out


def _predicted_shifts(molecule, nucleus) -> dict:
    """{atom id: [(shift, predicted ppm, label)]} in ascending ppm; 1H shifts are keyed by their carbon."""
    spectra = [s for s in molecule.spectra if s.nucleus == nucleus]
    out = {}
    if not spectra:
        return out
    prefix = "H" if nucleus == H1 else "C"
    for shift in max(spectra, key=lambda s: len(s.shifts)).shifts:
        if shift.ignored_auto or shift.ignored_user:
            continue
        value, _ = shift.get_value_error_pair(bjason.Molecule.CalcMethod.Best)
        value = None if value is None else float(value)
        mark = f"({shift.mark})" if shift.mark else ""
        for n in map(int, shift.nums):
            out.setdefault(n, []).append((shift, value, f"{prefix}{n + 1}{mark}"))
    for lst in out.values():
        lst.sort(key=lambda t: (t[1] is None, t[1]))
    return out


def _carbon_signals(c, h_ppms, c13_peaks, hc_peaks, eps, other_cs=(), used_c13_1d=True):
    """[(nmr_data index, peak, dim)] carrying this 13C value, a source description, and a note."""
    if c is None:
        return [], "", ""
    same_c = lambda p, dim: abs(float(p.pos[dim]) - c) <= eps
    hits = [(i, p, None) for i, p in c13_peaks if same_c(p, 0)]
    if hits:
        return hits[:1], f"13C 1D nmr_data[{hits[0][0]}]", ""
    direct = [(i, p, cd) for i, p, hd, cd in hc_peaks
              if same_c(p, cd) and any(abs(float(p.pos[hd]) - h) <= eps for h in h_ppms)]
    if direct:
        return direct[:1], f"2D direct nmr_data[{direct[0][0]}]", ""
    if not ASSIGN_LONG_RANGE_CARBONS or used_c13_1d:
        # If simpleNMR used a 1D 13C spectrum, its carbon values come from it,
        # so a miss there is a real problem, not something to bridge with HMBC peaks.
        return [], "", ""

    def belongs(p, cd):
        x = float(p.pos[cd])
        d = abs(x - c)
        return d <= LONG_RANGE_TOL and all(d < abs(x - o) for o in other_cs)

    long_range = [(i, p, cd) for i, p, hd, cd in hc_peaks if belongs(p, cd)]
    if not long_range:
        return [], "", ""
    mean = sum(float(p.pos[cd]) for _, p, cd in long_range) / len(long_range)
    if abs(mean - c) <= eps:
        note = f"via long-range peak (mean of {len(long_range)} = simpleNMR value)"
    else:
        note = f"via long-range peak (mean of {len(long_range)} is {mean:.4f}, simpleNMR {c:.4f})"
    # One shift -> one signal: use the peak closest to the value; JASON's merged
    # signals connect the carbon's other HMBC peaks to it.
    best = min(long_range, key=lambda t: abs(float(t[1].pos[t[2]]) - c))
    return [best], f"2D long-range nmr_data[{best[0]}]", note


def simplenmr_used_c13_1d(snmr: dict) -> bool:
    """True if the export says simpleNMR took carbon values from a 1D 13C spectrum."""
    for key, block in (snmr.get("oldjsondata") or {}).items():
        if isinstance(block, dict) and block.get("datatype") == "nmrspectrum":
            if block.get("experimenttype") == "C13_1D" or str(key).startswith("C13_1D"):
                return True
    return False


def build_assignments(document, snmr: dict, molecule_index: int = 0, eps: float = PPM_EPS):
    """Return (molecule, assignment tuples, AssignmentResult) for an open document.

    The tuples hold objects tied to the open file, so they must be used
    before the document is closed.
    """
    if snmr.get("peak_map"):
        return _build_by_id(document, snmr, molecule_index)
    nodes, elements, bonds = parse_simplenmr_export(snmr)
    molecule = document.mol_data[molecule_index]
    _check_atoms(molecule, elements, bonds)

    c13_peaks = _all_peaks(document, 1, ("13C",))
    hc_peaks = []                                   # (index, peak, h_dim, c_dim)
    for i, p in _all_peaks(document, 2, ("1H", "13C")):
        nuclides = tuple(document.nmr_data[i].spec_info.nuclides)
        hc_peaks.append((i, p, nuclides.index("1H"), nuclides.index("13C")))

    used_c13_1d = simplenmr_used_c13_1d(snmr)
    pred_c = _predicted_shifts(molecule, C13)
    pred_h = _predicted_shifts(molecule, H1)
    tuples, rows, seen = [], [], set()

    def add(tup):
        key = tuple(o if isinstance(o, int) or o is None else id(o) for o in tup)
        if key not in seen:            # a shift covering two symmetric atoms is listed under both
            seen.add(key)
            tuples.append(tup)

    def shifts_for(pred, atom_id):
        return [t for t in pred.get(atom_id, [])]

    for node in nodes:
        label, c, h_ppms = node["label"], node["c_ppm"], node["h_ppms"]

        # ---- 13C
        others = [n["c_ppm"] for n in nodes if n is not node and n["c_ppm"] is not None]
        signals, source, note = _carbon_signals(c, h_ppms, c13_peaks, hc_peaks, eps, others, used_c13_1d)
        c_shift_lists = [shifts_for(pred_c, a) for a in node["atom_ids"]]
        if c is None:
            rows.append(ReportRow(label, "13C", "", None, None, "", "no 13C value in simpleNMR"))
        elif not signals:
            rows.append(ReportRow(label, "13C", "", None, c, "", "no exact 13C peak"))
        elif not any(c_shift_lists):
            rows.append(ReportRow(label, "13C", "", None, c, source, "no predicted 13C shift"))
        else:
            status = "assigned" + (f" {note}" if note else "")
            done = set()
            for shift, value, shift_label in (t for lst in c_shift_lists for t in lst):
                if id(shift) in done:
                    continue
                done.add(id(shift))
                for i, peak, dim in signals:
                    add((shift, peak) if dim is None else (shift, peak, dim))
                rows.append(ReportRow(label, "13C", shift_label, value, c, source, status))

        # ---- 1H: exact match on the 1H axis; closest 13C coordinate picks the HSQC peak
        h_signals = []
        for h in h_ppms:
            hits = [t for t in hc_peaks if abs(float(t[1].pos[t[2]]) - h) <= eps]
            if not hits:
                rows.append(ReportRow(label, "1H", "", None, h, "", "no exact 1H peak"))
                continue
            if c is not None:
                hits.sort(key=lambda t: abs(float(t[1].pos[t[3]]) - c))
            h_signals.append((h, *hits[0]))
        if not h_signals:
            continue
        h_shift_lists = [shifts_for(pred_h, a) for a in node["atom_ids"]]
        if not any(h_shift_lists):
            rows.append(ReportRow(label, "1H", "", None, h_signals[0][0], "", "no predicted 1H shift"))
            continue
        done = set()
        for all_h_shifts in h_shift_lists:            # pair per atom, so symmetric CH2s pair correctly
            h_shifts = [t for t in all_h_shifts if id(t[0]) not in done]
            done.update(id(t[0]) for t in h_shifts)
            sig = h_signals * len(h_shifts) if len(h_signals) == 1 and len(h_shifts) > 1 else h_signals
            for (shift, value, shift_label), (h, i, peak, h_dim, _c_dim) in zip(h_shifts, sig):
                add((shift, peak, h_dim))
                rows.append(ReportRow(label, "1H", shift_label, value, h, f"nmr_data[{i}]", "assigned"))

    _check_one_signal_per_shift(tuples)

    n_sym = sum(len(n["atom_ids"]) - 1 for n in nodes)
    return molecule, tuples, AssignmentResult(len(tuples), len(nodes) + n_sym, rows, n_sym)


def _peaks_by_id(document) -> dict:
    """{JEOL ID: (peak, nuclides of its spectrum)} for every peak in the document."""
    out = {}
    for spectrum in document.nmr_data:
        nuclides = tuple(spectrum.spec_info.nuclides)
        try:
            for peak in spectrum.peaks:
                out[str(peak.id)] = (peak, nuclides)
        except KeyError:        # spectrum without a peak table
            pass
    return out


def _build_by_id(document, snmr: dict, molecule_index: int = 0):
    """Route 1: assignments and peak moves from the export's peak_map (see module docstring)."""
    nodes, elements, bonds = parse_simplenmr_export(snmr)
    molecule = document.mol_data[molecule_index]
    _check_atoms(molecule, elements, bonds)

    peak_map = snmr["peak_map"]
    peaks = _peaks_by_id(document)
    pred_c = _predicted_shifts(molecule, C13)
    pred_h = _predicted_shifts(molecule, H1)
    tuples, rows, seen = [], [], set()

    def add(tup):
        key = tuple(o if isinstance(o, int) or o is None else id(o) for o in tup)
        if key not in seen:
            seen.add(key)
            tuples.append(tup)

    def lookup(jason_id, nucleus):
        """(peak, axis of `nucleus` or None for a 1D peak), or None if the ID isn't in the document."""
        if jason_id not in peaks:
            return None
        peak, nuclides = peaks[jason_id]
        return peak, (None if len(nuclides) == 1 else nuclides.index(nucleus))

    for node in nodes:
        label = node["label"]
        group = peak_map["nodes"].get(node["node_id"])
        if group is None:
            rows.append(ReportRow(label, "13C", "", None, node["c_ppm"], "", "not in peak map"))
            continue

        # ---- carbon: 1D 13C peak, else HSQC 13C axis, else an HMBC peak's 13C axis
        if group.get("c13_1d"):
            c_id, source = group["c13_1d"], "13C 1D"
        elif group.get("hsqc"):
            c_id, source = group["hsqc"][0]["jason_id"], "HSQC"
        else:
            hmbc = [c for c in peak_map["correlations"]
                    if c["kind"] == "HMBC" and c["delta1_node"] == node["node_id"]]
            c_id, source = (hmbc[0]["jason_id"], "HMBC") if hmbc else (None, "")
        c_signal = lookup(c_id, "13C") if c_id else None
        c_shift_lists = [pred_c.get(a, []) for a in node["atom_ids"]]
        if c_signal is None:
            rows.append(ReportRow(label, "13C", "", None, node["c_ppm"], source,
                                  f"JEOL ID {c_id} not found in the document" if c_id else "no carbon peak in peak map"))
        elif not any(c_shift_lists):
            rows.append(ReportRow(label, "13C", "", None, node["c_ppm"], source, "no predicted 13C shift"))
        else:
            peak, dim = c_signal
            for shift, value, shift_label in (t for lst in c_shift_lists for t in lst):
                add((shift, peak) if dim is None else (shift, peak, dim))
                rows.append(ReportRow(label, "13C", shift_label, value, node["c_ppm"], source, "assigned"))

        # ---- protons: HSQC peak(s) on the 1H axis; per atom, so symmetric CH2s pair correctly
        hsqc = sorted(group.get("hsqc", []), key=lambda h: h["h_ppm"])
        h_signals = []
        for h in hsqc:
            signal = lookup(h["jason_id"], "1H")
            if signal is None:
                rows.append(ReportRow(label, "1H", "", None, h["h_ppm"], "HSQC",
                                      f"JEOL ID {h['jason_id']} not found in the document"))
            else:
                h_signals.append((h["h_ppm"], *signal))
        if not h_signals:
            continue
        for atom in node["atom_ids"]:
            h_shifts = pred_h.get(atom, [])
            if not h_shifts:
                rows.append(ReportRow(label, "1H", "", None, h_signals[0][0], "HSQC", "no predicted 1H shift"))
                continue
            # one HSQC peak: every proton shift on it; two (diastereotopic CH2): pair in ascending ppm
            signals = h_signals * len(h_shifts) if len(h_signals) == 1 else h_signals
            for (shift, value, shift_label), (h_ppm, peak, dim) in zip(h_shifts, signals):
                add((shift, peak, dim))
                rows.append(ReportRow(label, "1H", shift_label, value, h_ppm, "HSQC", "assigned"))

    # ---- HMBC/COSY peaks: move to simpleNMR's snapped positions (delta2 = Pos[0], delta1 = Pos[1])
    moves = []
    for c in peak_map.get("correlations", []):
        if c["jason_id"] not in peaks:
            rows.append(ReportRow(c["kind"], c["kind"], "", None, None, "",
                                  f"JEOL ID {c['jason_id']} not found in the document"))
            continue
        peak = peaks[c["jason_id"]][0]
        original = list(peak.pos)
        new = list(original)
        new[0], new[1] = c["delta2"], c["delta1"]
        moves.append(PeakMove(c["kind"], c["jason_id"], peak, original, new))

    _check_one_signal_per_shift(tuples)
    n_sym = sum(len(n["atom_ids"]) - 1 for n in nodes)
    result = AssignmentResult(len(tuples), len(nodes) + n_sym, rows, n_sym, moves)
    return molecule, tuples, result


def _check_one_signal_per_shift(tuples):
    """JASON rejects a shift assigned to more than one signal; fail early with a clear message."""
    seen = {}
    for t in tuples:
        key = id(t[0])
        if key in seen and seen[key] != (id(t[1]), t[2] if len(t) > 2 else None):
            raise ValueError("internal error: a predicted shift was given more than one signal")
        seen[key] = (id(t[1]), t[2] if len(t) > 2 else None)


# ---------------------------------------------------------------- file-level entry points

def preview_assignments(jjh5_path: Path, snmr: dict, **kwargs) -> AssignmentResult:
    """Build the assignments read-only and report what would be written."""
    with bjason.Document(str(jjh5_path), mode="r") as doc:
        _, _, result = build_assignments(doc, snmr, **kwargs)
    return result


def write_assignments(jjh5_path: Path, snmr: dict, backup: bool = False, **kwargs) -> AssignmentResult:
    """Write the assignments into jjh5_path in place.

    With backup=True the file is first copied to
    ``<stem>.<timestamp>.bak.jjh5`` alongside it. On the ID route, HMBC/COSY
    peaks are moved and their original positions written to
    ``<stem>_peak_moves.json`` alongside it.
    """
    if not hasattr(bjason.Molecule, "set_assignments"):
        raise RuntimeError(
            "This beautifuljason has no Molecule.set_assignments - "
            "install beautifuljason >= 1.3.0b1."
        )
    jjh5_path = Path(jjh5_path)
    if backup:
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        shutil.copy2(jjh5_path, jjh5_path.with_name(f"{jjh5_path.stem}.{stamp}.bak{jjh5_path.suffix}"))
    # Rebuild inside the r+ session: shift and peak objects belong to the open file.
    with bjason.Document(str(jjh5_path), mode="r+") as doc:
        molecule, tuples, result = build_assignments(doc, snmr, **kwargs)
        # Move the HMBC/COSY peaks first; JASON regroups the signals when it loads the file.
        for move in result.moves:
            move.peak.h5_group.attrs["Pos"] = move.new
        molecule.set_assignments(tuples)
    if result.moves:
        # Record the edit to the peak table so it is visible and reversible.
        log = jjh5_path.with_name(f"{jjh5_path.stem}_peak_moves.json")
        log.write_text(json.dumps([m.log_entry() for m in result.moves], indent=2))
    return result
