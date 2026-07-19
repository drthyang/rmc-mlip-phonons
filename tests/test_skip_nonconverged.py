"""Unit tests for read_moves_generated / drop_nonconverged.

RMCProfile writes ``Number of moves generated: N`` in the rmc6f header; a
config with N = 0 is a run that never started (non-converged). Files without
the line (synthetic fixtures, other generators) must be kept.
"""

from __future__ import annotations

from conftest import write_rmc6f

HEADER = """
    (Version 6f format configuration file)
    Number of moves generated:           {moves}
    Number of atoms: 1
    Supercell dimensions: 1 1 1
    Cell (Ang/deg): 4.0 4.0 4.0 90 90 90
    Atoms:
        1  Cu  0.0 0.0 0.0
"""

NO_MOVES_HEADER = """
    (Version 6 format configuration file)
    Number of atoms: 1
    Supercell dimensions: 1 1 1
    Cell (Ang/deg): 4.0 4.0 4.0 90 90 90
    Atoms:
        1  Cu  0.0 0.0 0.0
"""


def test_read_moves_generated_parses_count(tmp_path, m1):
    f = write_rmc6f(tmp_path / "a.rmc6f", HEADER.format(moves=2230804))
    assert m1.read_moves_generated(f) == 2230804


def test_read_moves_generated_zero(tmp_path, m1):
    f = write_rmc6f(tmp_path / "z.rmc6f", HEADER.format(moves=0))
    assert m1.read_moves_generated(f) == 0


def test_read_moves_generated_absent_returns_none(tmp_path, m1):
    f = write_rmc6f(tmp_path / "n.rmc6f", NO_MOVES_HEADER)
    assert m1.read_moves_generated(f) is None


def test_read_moves_stops_at_atoms_section(tmp_path, m1):
    """A number in the atoms section must not be misread as a move count."""
    text = NO_MOVES_HEADER + "    2  Cu  0.5 0.5 0.5\n"
    f = write_rmc6f(tmp_path / "s.rmc6f", text)
    assert m1.read_moves_generated(f) is None


def test_drop_nonconverged_partitions_and_preserves_order(tmp_path, m1):
    good1 = write_rmc6f(tmp_path / "c1.rmc6f", HEADER.format(moves=100))
    bad = write_rmc6f(tmp_path / "c2.rmc6f", HEADER.format(moves=0))
    good2 = write_rmc6f(tmp_path / "c3.rmc6f", HEADER.format(moves=200))
    noheader = write_rmc6f(tmp_path / "c4.rmc6f", NO_MOVES_HEADER)

    kept, dropped = m1.drop_nonconverged([good1, bad, good2, noheader])
    assert kept == [good1, good2, noheader]  # absent header -> kept
    assert dropped == [bad]


def test_drop_nonconverged_keeps_synthetic_fixture(tmp_path, m1):
    """The Cu fcc fixture generator writes no moves line -> nothing dropped."""
    from fixtures.make_synthetic_ensemble import make_fcc_cu_ensemble

    paths = make_fcc_cu_ensemble(tmp_path, n_configs=3, seed=0)
    kept, dropped = m1.drop_nonconverged(paths)
    assert kept == paths
    assert dropped == []
