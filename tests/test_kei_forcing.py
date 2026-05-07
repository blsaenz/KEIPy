"""Tests for `kei_forcing` dataset assembly (pure Python + xarray)."""

from __future__ import annotations

import numpy as np

import kei as keipy

from .helpers_forcing import minimal_forcing_dict


def test_kei_forcing_builds_dataset_from_dict():
    f_dict = minimal_forcing_dict(nz=48, nt=5)
    ds = keipy.kei_forcing(f_dict=f_dict)

    assert "f_time" in ds.variables
    assert "zm" in ds.variables
    assert "dm" in ds.variables
    assert "hm" in ds.variables
    assert ds.sizes["zm"] == 48
    assert ds.sizes["f_time"] == 5

    for name in keipy.forcing_idx:
        assert name in ds.variables
        assert ds[name].dims == ("f_time",)

    for v in keipy.init_vars_ocn + keipy.init_vars_eco:
        assert v in ds.variables
        assert ds[v].dims == ("zm",)


def test_kei_runtime_yaml_has_driver_keys_under_kei_common():
    """Driver scalars live under ``kei_common`` in YAML, matching ``kei.kei_common``."""
    from pathlib import Path

    yaml_path = Path(__file__).resolve().parent.parent / "kei_runtime_params.yml"
    doc = keipy._load_runtime_yaml(str(yaml_path))
    kc = doc["kei_common"]
    for k in ("dtsec", "dlon", "dlat", "lice", "leco", "lsw"):
        assert k in kc
