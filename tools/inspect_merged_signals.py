"""
Read-only diagnostic: JASON's merged signals, which peaks each contains,
which ones receive an assignment, and which correlations link them.

Usage (py313arm environment, beautifuljason >= 1.3.0b1):

    python tools/inspect_merged_signals.py <input .jjh5 or .bak.jjh5> <html/input_py.json>
"""

import json
import sys

import beautifuljason as bjason

from simplenmrjeol.assignment_writer import build_assignments

H1 = bjason.Molecule.Atom.NuclType.H1
C13 = bjason.Molecule.Atom.NuclType.C13


def public(obj):
    return [a for a in dir(obj) if not a.startswith("_")]


def peaks_of(spectrum):
    out = []
    for kind in ("peaks", "multiplets"):
        try:
            out += [(kind, p) for p in getattr(spectrum, kind)]
        except (KeyError, AttributeError):
            pass
    return out


def fmt_pos(p):
    return "(" + ", ".join(f"{float(x):.4f}" for x in p.pos) + ")"


def main(jjh5, export):
    snmr = json.loads(open(export, encoding="utf-8").read())
    with bjason.Document(jjh5, mode="r") as doc:
        molecule, tuples, result = build_assignments(doc, snmr)
        print(result.summary_text(), "\n")

        print("Spectra")
        for i, s in enumerate(doc.nmr_data):
            print(f"  nmr_data[{i}] {s.ndim}D {tuple(s.spec_info.nuclides)}  {len(peaks_of(s))} peaks/multiplets")

        merged_1d = list(molecule.merged_1d_signals)
        merged_2d = list(molecule.merged_2d_signals)
        print(f"\n{len(merged_1d)} merged 1D signals, {len(merged_2d)} merged 2D signals")
        if merged_1d:
            print("  merged 1D attributes:", public(merged_1d[0]))
        if merged_2d:
            print("  merged 2D attributes:", public(merged_2d[0]))

        # what the writer assigns: (shift label, nmr_data index, dim, peak)
        index_of = {}
        for i, s in enumerate(doc.nmr_data):
            for _, p in peaks_of(s):
                index_of[id(p)] = i
        assigned = []
        for t in tuples:
            shift, sig = t[0], t[1]
            nuc = "H" if shift in [x for sp in molecule.spectra if sp.nucleus == H1 for x in sp.shifts] else "C"
            label = nuc + "/".join(str(int(n) + 1) for n in shift.nums) + (f"({shift.mark})" if shift.mark else "")
            assigned.append((label, sig, t[2] if len(t) > 2 else None))

        short = {}
        print("\nMerged 1D signals: contents and assignments")
        for k, m in enumerate(sorted(merged_1d, key=lambda m: (m.nucleus.name, m.position))):
            short[m.id] = f"M{k}"
            got = [f"{lab}@nmr_data[{index_of.get(id(sig), '?')}]{'' if d is None else f'dim{d}'}"
                   for lab, sig, d in assigned if m.references(sig)]
            print(f"  M{k:<3} {m.nucleus.name:4} {m.position:9.4f}  assigned: {', '.join(got) or 'NONE'}")
            for i, s in enumerate(doc.nmr_data):
                for kind, p in peaks_of(s):
                    try:
                        if m.references(p):
                            print(f"         contains nmr_data[{i}] {kind[:-1]} {fmt_pos(p)}")
                    except Exception as e:      # noqa: BLE001 - diagnostic only
                        print(f"         references() failed on nmr_data[{i}]: {e!r}")
                        break

        print("\nMerged 2D signals")
        for m2 in merged_2d:
            ends = []
            for dim in range(2):
                parts = []
                for nuc in (H1, C13):
                    for sid in m2.linked_signal_ids(dim, nuc):
                        parts.append(f"{nuc.name}:{short.get(sid, str(sid))}")
                ends.append(" ".join(parts) or "-")
            extra = f" bonds {m2.min_distance}-{m2.max_distance}" if getattr(m2, "min_distance", None) is not None else ""
            extra += " IGNORED" if getattr(m2, "ignored", False) else ""
            pos = getattr(m2, "position", None) or getattr(m2, "pos", None)
            pos = f" at {tuple(round(float(x), 4) for x in pos)}" if pos is not None else ""
            print(f"  [{ends[0]}]  <->  [{ends[1]}]{extra}{pos}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
