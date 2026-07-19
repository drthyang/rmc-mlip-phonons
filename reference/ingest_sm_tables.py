#!/usr/bin/env python3
"""Ingest the GaTa4Se8 supplemental tables into gts_mode_patterns.json.

Source: Supplemental Material of Yang et al., Phys. Rev. Research 4, 033123
(2022) — local file (git-ignored) at data/GaTa4Se8_SM.pdf. Extracted tables:

    II    crystallography-refined P-4̄2₁m structure (20 sites, SHELX labels)
    III   PDF-refined local structure (same labeling as II)
    IV    total normalized atomic displacements at 10 K (AMPLIMODES labels —
          note: a DIFFERENT site labeling than Table II)
    V–X   per-irrep normalized displacement patterns (Γ₁, Γ₃, Δ₃*, X₃, X₅,
          W₄), sharing Table IV's labeling
          (*the SM labels the Δ mode "Δ₃" while main-text Table II calls it
          "Δ₁" — recorded verbatim with this note)

Units: fractional coordinates of the tetragonal cell (a = 10.34370,
c = 20.6878 Å); |u| and mode amplitudes in Å. Published mode amplitudes
(main-text Table II, normalized w.r.t. the primitive F-4̄3m cell):
Γ₁ 0.0061, Γ₃ 0.0020, Δ 0.0212, X₃ 0.0719, X₅ 0.1196, W₄ 0.0260 Å.

The empirical normalization of Tables V–X is *not assumed*: the test suite
decomposes Table IV onto the patterns and checks the coefficient ratios
against the published amplitudes.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PDF = REPO / "data/GaTa4Se8_SM.pdf"
OUT = Path(__file__).resolve().parent / "gts_mode_patterns.json"

CELL = {"a_A": 10.34370, "c_A": 20.6878, "spacegroup": "P-42_1m (113)"}
PUBLISHED_AMPLITUDES_A = {  # main-text Table II of PRR 4, 033123
    "G1": 0.0061, "G3": 0.0020, "D": 0.0212,
    "X3": 0.0719, "X5": 0.1196, "W4": 0.0260,
}
MODE_TABLES = {  # SM table -> mode key, irrep label as printed, k-vector
    "V": ("G1", "Gamma1", [0, 0, 0]),
    "VI": ("G3", "Gamma3", [0, 0, 0]),
    "VII": ("D", "Delta3 (main text: Delta1)", [0, 0, 0.5]),
    "VIII": ("X3", "X3", [0, 1, 0]),
    "IX": ("X5", "X5", [0, 1, 0]),
    "X": ("W4", "W4", [0.5, 1, 0]),
}

ROW = re.compile(
    r"^(?P<label>(?:Ta|Se|Ga)\d+)\s+(?P<wp>4e|8f)\s+(?P<rest>[-0-9. ()]+)$")

# Documented corrections to the *printed* SM (applied at parse time, recorded
# in the output JSON; flagged for author verification against AMPLIMODES):
CORRECTIONS = [
    {
        "table": "VIII", "label": "Se3", "component": 2,
        "printed": -0.0053, "corrected": 0.0053,
        "justification": "Least-squares reconstruction of Table IV from the "
            "six mode patterns leaves a single outlier at Se3 (residual > |u|)"
            " that a +z sign resolves; and in every other mode table Se3/Se6 "
            "are antisymmetric partners, while Table VIII prints them "
            "identical. Printed -0.0053 is a probable sign typo.",
    },
]


def _numbers(s: str):
    """Floats from a row remainder, dropping parenthesized esds."""
    return [float(t) for t in re.findall(r"-?\d+\.\d+", re.sub(r"\(\d+\)", "", s))]


def parse(pdf_path: Path):
    from pypdf import PdfReader

    text = "\n".join(p.extract_text() for p in PdfReader(str(pdf_path)).pages)
    # split into table chunks: "TABLE <ROMAN>." headers
    chunks = re.split(r"TABLE\s+([IVX]+)\.", text)
    tables = {}
    for i in range(1, len(chunks) - 1, 2):
        roman = chunks[i]
        rows = {}
        for line in chunks[i + 1].splitlines():
            m = ROW.match(line.strip())
            if m:
                nums = _numbers(m.group("rest"))
                rows[m.group("label")] = {"wp": m.group("wp"), "values": nums}
        tables[roman] = rows
    return tables


def main():
    if not PDF.is_file():
        sys.exit(f"missing {PDF} — copy the supplemental PDF there first")
    tables = parse(PDF)

    for corr in CORRECTIONS:
        row = tables[corr["table"]][corr["label"]]
        assert row["values"][corr["component"]] == corr["printed"], corr
        row["values"][corr["component"]] = corr["corrected"]
        print(f"applied correction: table {corr['table']} {corr['label']} "
              f"component {corr['component']}: {corr['printed']} -> "
              f"{corr['corrected']}")

    def structure(roman):
        return {lab: {"wp": r["wp"], "xyz": r["values"][:3]}
                for lab, r in tables[roman].items()}

    out = {
        "source": {
            "citation": "Yang et al., Phys. Rev. Research 4, 033123 (2022), "
                        "Supplemental Material",
            "file": "data/GaTa4Se8_SM.pdf",
            "note": "Tables IV-X use AMPLIMODES site labels, which differ "
                    "from Table II/III SHELX labels; the mapping is resolved "
                    "numerically in mode_project.py, never by name.",
        },
        "cell": CELL,
        "corrections_applied": CORRECTIONS,
        "published_amplitudes_A": PUBLISHED_AMPLITUDES_A,
        "table_II_refined": structure("II"),
        "table_III_pdf_refined": structure("III"),
        "table_IV_total": {
            lab: {"wp": r["wp"], "u_frac": r["values"][:3],
                  "u_abs_A": r["values"][3]}
            for lab, r in tables["IV"].items()},
        "modes": {},
    }
    for roman, (key, irrep, k) in MODE_TABLES.items():
        out["modes"][key] = {
            "sm_table": roman, "irrep": irrep, "k_conventional_cubic": k,
            "amplitude_published_A": PUBLISHED_AMPLITUDES_A[key],
            "pattern_frac": {lab: {"wp": r["wp"], "d_frac": r["values"][:3]}
                             for lab, r in tables[roman].items()},
        }

    counts = {roman: len(tables[roman]) for roman in tables}
    print("rows parsed per table:", counts)
    OUT.write_text(json.dumps(out, indent=1))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
