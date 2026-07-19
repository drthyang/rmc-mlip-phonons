"""Unit tests for milestone1_bands.select_configs (ensemble subsampling).

--stride decimates in order; --max-configs then caps via a seeded random
sample; the whole thing is deterministic in --seed and order-preserving.
select_configs never touches the filesystem, so plain Path names suffice.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _files(n):
    return [Path(f"c{i:03d}.rmc6f") for i in range(n)]


def test_no_subsampling_returns_all(m1):
    fs = _files(5)
    assert m1.select_configs(fs) == fs


def test_stride_decimates_in_order(m1):
    fs = _files(10)
    assert m1.select_configs(fs, stride=3) == fs[::3]  # 0,3,6,9


def test_max_configs_caps_count_without_duplicates(m1):
    fs = _files(100)
    sel = m1.select_configs(fs, max_configs=10, seed=0)
    assert len(sel) == 10
    assert len(set(sel)) == 10          # replace=False
    assert set(sel) <= set(fs)          # subset of the inputs


def test_max_configs_preserves_input_order(m1):
    fs = _files(100)
    sel = m1.select_configs(fs, max_configs=20, seed=1)
    idx = [fs.index(p) for p in sel]
    assert idx == sorted(idx)


def test_same_seed_is_deterministic(m1):
    fs = _files(100)
    assert (m1.select_configs(fs, max_configs=15, seed=42)
            == m1.select_configs(fs, max_configs=15, seed=42))


def test_different_seed_changes_selection(m1):
    fs = _files(200)
    assert (m1.select_configs(fs, max_configs=20, seed=1)
            != m1.select_configs(fs, max_configs=20, seed=2))


def test_stride_then_max_configs_compose(m1):
    fs = _files(100)
    sel = m1.select_configs(fs, stride=2, max_configs=10, seed=0)
    assert len(sel) == 10
    # everything sampled came from the strided (even-index) subset
    assert all(fs.index(p) % 2 == 0 for p in sel)


def test_cap_not_binding_returns_all_after_stride(m1):
    fs = _files(5)
    assert m1.select_configs(fs, max_configs=10) == fs      # cap > len
    assert m1.select_configs(fs, max_configs=5) == fs       # cap == len
    assert m1.select_configs(fs, stride=2, max_configs=10) == fs[::2]


@pytest.mark.parametrize("kwargs, msg", [
    ({"stride": 0}, "stride"),
    ({"stride": -1}, "stride"),
    ({"max_configs": 0}, "max-configs"),
])
def test_invalid_arguments_raise(m1, kwargs, msg):
    with pytest.raises(SystemExit, match=msg):
        m1.select_configs(_files(4), **kwargs)
