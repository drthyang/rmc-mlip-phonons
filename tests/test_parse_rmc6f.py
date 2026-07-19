"""Unit tests for milestone1_bands.parse_rmc6f (tolerant .rmc6f reader).

Covers the format variations called out in CLAUDE.md immediate task 1:
with/without site-id columns, bracketed labels, and
``Cell (Ang/deg)`` vs ``Lattice vectors`` cell headers.
"""

from __future__ import annotations

import numpy as np
import pytest

from conftest import write_rmc6f

# A supercell with explicit lattice vectors, bracketed labels, and a trailing
# site-id column. Two unit cells along x (dims 2 1 1) so folding is testable.
RMC6F_WITH_SITE_IDS = """
    (Version 6 format configuration file)
    Metadata date: 2026-07-19
    Number of atoms: 2
    Supercell dimensions:        2 1 1
    Lattice vectors (Ang):
        7.20  0.00  0.00
        0.00  3.60  0.00
        0.00  0.00  3.60
    Atoms:
         1   Cu   [1]    0.10  0.20  0.30    1  0  0  0
         2   Cu   [1]    0.60  0.20  0.30    1  1  0  0
"""

# Same structure via the Cell (Ang/deg) header, no brackets, no site ids,
# and one deliberately negative coordinate to exercise the % 1.0 wrap.
RMC6F_CELL_HEADER_NO_IDS = """
    (Version 6 format configuration file)
    Number of atoms: 2
    Supercell dimensions:        2 1 1
    Cell (Ang/deg):   7.20 3.60 3.60 90.0 90.0 90.0
    Atoms:
         1   Cu   0.10   0.20   0.30
         2   Cu   0.60  -0.80   0.30
"""


def test_parse_with_site_ids(tmp_path, m1):
    cfg = m1.parse_rmc6f(write_rmc6f(tmp_path / "a.rmc6f", RMC6F_WITH_SITE_IDS))

    assert cfg["elements"] == ["Cu", "Cu"]
    np.testing.assert_array_equal(cfg["dims"], [2, 1, 1])
    # Explicit lattice-vector rows read verbatim.
    np.testing.assert_allclose(
        cfg["cell"], [[7.2, 0, 0], [0, 3.6, 0], [0, 0, 3.6]], atol=1e-9
    )
    np.testing.assert_allclose(cfg["frac"], [[0.10, 0.20, 0.30], [0.60, 0.20, 0.30]])
    # Site-id column present -> the first trailing integer is the site id.
    assert cfg["site_ids"] is not None
    np.testing.assert_array_equal(cfg["site_ids"], [1, 1])


def test_parse_without_site_ids_and_cell_header(tmp_path, m1):
    cfg = m1.parse_rmc6f(
        write_rmc6f(tmp_path / "b.rmc6f", RMC6F_CELL_HEADER_NO_IDS)
    )

    assert cfg["site_ids"] is None  # no trailing integer column
    np.testing.assert_array_equal(cfg["dims"], [2, 1, 1])
    # Cell built from a b c alpha beta gamma matches the lattice-vector form.
    np.testing.assert_allclose(
        cfg["cell"], [[7.2, 0, 0], [0, 3.6, 0], [0, 0, 3.6]], atol=1e-6
    )
    # -0.80 wraps into [0, 1): -0.80 % 1 == 0.20.
    np.testing.assert_allclose(cfg["frac"], [[0.10, 0.20, 0.30], [0.60, 0.20, 0.30]])


def test_cell_header_and_lattice_vectors_agree(tmp_path, m1):
    """The two supported cell headers must yield the same lattice."""
    a = m1.parse_rmc6f(write_rmc6f(tmp_path / "a.rmc6f", RMC6F_WITH_SITE_IDS))
    b = m1.parse_rmc6f(write_rmc6f(tmp_path / "b.rmc6f", RMC6F_CELL_HEADER_NO_IDS))
    np.testing.assert_allclose(a["cell"], b["cell"], atol=1e-6)


def test_bracketed_label_is_skipped(tmp_path, m1):
    """A ``[label]`` token between element and coordinates is not a coordinate."""
    text = """
        Supercell dimensions: 1 1 1
        Lattice vectors (Ang):
            3.60 0.00 0.00
            0.00 3.60 0.00
            0.00 0.00 3.60
        Atoms:
            1  Cu  [1]  0.11  0.22  0.33
    """
    cfg = m1.parse_rmc6f(write_rmc6f(tmp_path / "brk.rmc6f", text))
    assert cfg["elements"] == ["Cu"]
    np.testing.assert_allclose(cfg["frac"], [[0.11, 0.22, 0.33]])
    assert cfg["site_ids"] is None


def test_element_cleaning_from_ionic_labels(tmp_path, m1):
    """Charged/upper-case element tokens reduce to a proper symbol."""
    text = """
        Supercell dimensions: 1 1 1
        Cell (Ang/deg): 4.0 4.0 4.0 90 90 90
        Atoms:
            1  O2-  0.00 0.00 0.00
            2  MN   0.50 0.50 0.50
    """
    cfg = m1.parse_rmc6f(write_rmc6f(tmp_path / "ion.rmc6f", text))
    assert cfg["elements"] == ["O", "Mn"]


@pytest.mark.parametrize(
    "text, msg",
    [
        (
            "Supercell dimensions: 1 1 1\nAtoms:\n  1 Cu 0.0 0.0 0.0\n",
            "no lattice information",
        ),
        (
            "Cell (Ang/deg): 4 4 4 90 90 90\nAtoms:\n  1 Cu 0.0 0.0 0.0\n",
            "Supercell",
        ),
        (
            "Supercell dimensions: 1 1 1\nCell (Ang/deg): 4 4 4 90 90 90\nAtoms:\n",
            "no atom lines",
        ),
    ],
)
def test_missing_sections_raise(tmp_path, m1, text, msg):
    with pytest.raises(ValueError, match=msg):
        m1.parse_rmc6f(write_rmc6f(tmp_path / "bad.rmc6f", text))
