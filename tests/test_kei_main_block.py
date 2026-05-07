"""
Integration test mirroring `kei.py` ``if __name__ == '__main__':`` (MACMODS / lsw=1).

Requires ``test_data/kf_200_100_2000.nc``, reference ``test_data/test_output.nc``,
and a built ``f90.kei`` extension. The run is lengthy (~768 hourly steps); use
``pytest -m slow`` to run only this module or ``pytest -m 'not slow'`` to skip it.

Cross-environment fidelity (why e.g. py311 can match ``test_output.nc`` but py313
does not):

- **Committed reference**: ``test_output.nc`` is a frozen run. Bit-level or tight
  float agreement across Python/OS/BLAS stacks is not guaranteed for a nonlinear
  month-long integration (`T` couples ice, KPP, tracers).

- **Fortran extension link**: Py311 vs py313 builds often link **different LAPACK**
  (e.g. macOS Accelerate vs conda-forge OpenBLAS / ``liblapack``). ``SGTSV`` and related
  solvers accumulate small differences each step (~10³ steps ⇒ order 0.1–1 °C).

- **Threading**: OpenBLAS / MKL can reorder work; set ``OPENBLAS_NUM_THREADS=1``,
  ``OMP_NUM_THREADS=1``, ``MKL_NUM_THREADS=1`` when debugging drift.

- **Python stack**: Newer pandas/xarray/cftime may shift **time resampling**
  (``compute`` ``resample(…).interpolate()`` → ``nt`` or boundary forcing). If
  ``nt`` or the last ``f_time`` label differs between envs, fix the upstream grid
  first.

To refresh the reference after changing toolchains: rerun once in your canonical env
and replace ``test_data/test_output.nc``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

_REPO = Path(__file__).resolve().parent.parent
_TEST_DATA = _REPO / "test_data"
_INPUT_NC = _TEST_DATA / "kf_200_100_2000.nc"
_REF_NC = _TEST_DATA / "test_output.nc"


def _env_digest_for_comparison() -> str:
    bits = [
        f"Python {sys.version.split()[0]}",
        f"numpy {np.__version__}",
    ]
    try:
        import pandas as pd  # noqa: PLC0415

        bits.append(f"pandas {pd.__version__}")
    except ImportError:
        pass
    try:
        import xarray as xr  # noqa: PLC0415

        bits.append(f"xarray {xr.__version__}")
    except ImportError:
        pass
    return ", ".join(bits)


def _run_main_block_workflow(out_base: Path) -> tuple[Path, Path]:
    """Mirror ``kei.py`` __main__: MACMODS forcing fields + short date range."""
    import kei as keipy

    if keipy.kei is None:
        pytest.skip(
            "Fortran extension missing — build with `make -C f90 kei` for this environment."
        )

    if not _INPUT_NC.is_file():
        pytest.skip(f"Input NetCDF not found: {_INPUT_NC}")
    if not _REF_NC.is_file():
        pytest.skip(f"Reference NetCDF not found: {_REF_NC}")

    kf_ds = keipy.kei_forcing(
        nc_file=str(_INPUT_NC),
        start_date="2000-01-01",
        freq="h",
        legacy_nc=True,
    )
    f_time_dim = kf_ds.sizes["f_time"]
    kf_ds["swh"] = ("f_time"), np.full(f_time_dim, 0.5)
    kf_ds["mwp"] = ("f_time"), np.full(f_time_dim, 30.0)
    kf_ds["cmag"] = ("f_time"), np.full(f_time_dim, 0.05)

    k = keipy.kei_simulation(
        kf_ds,
        t_start="2000-01-15",
        t_end="2000-02-15",
    )
    out_base.mkdir(parents=True, exist_ok=True)
    k.apply_runtime_yaml('kei_runtime_params.yml')
    k.compute(
        str(out_base),
        run_name="test_output",
        yaml_overrides={"kei_common": {"lsw": 1}},
    )

    run_dirs = sorted(out_base.glob("test_output_*"))
    assert len(run_dirs) == 1, f"expected one output directory, got {run_dirs!r}"
    run_dir = run_dirs[0]
    out_nc = run_dir / "test_output.nc"
    params_csv = run_dir / "test_output_keipy_params.csv"
    assert out_nc.is_file(), f"missing {out_nc}"
    assert params_csv.is_file(), f"missing {params_csv}"
    return out_nc, params_csv


@pytest.fixture(scope="module")
def main_block_paths(tmp_path_factory):
    base = tmp_path_factory.mktemp("kei_main_block")
    return _run_main_block_workflow(base)


@pytest.mark.fortran
@pytest.mark.slow
def test_main_block_writes_params_csv(main_block_paths):
    _, params_csv = main_block_paths
    text = params_csv.read_text(encoding="utf-8")
    assert "lsw" in text


def _array_last_time_close(
    var: str,
    ref_nc: Path,
    run_nc: Path,
    *,
    rtol: float = 1e-5,
    atol: float = 1e-4,
) -> None:
    with xr.open_dataset(ref_nc) as ref, xr.open_dataset(run_nc) as run:
        assert var in ref.data_vars, f"{var} not in reference"
        assert var in run.data_vars, f"{var} not in run output"
        # f_time length can differ (forcing resample); compare the same clock time as the ref end.
        t_end = ref["f_time"].isel(f_time=-1)
        da_r = ref[var].sel(f_time=t_end)
        da_y = run[var].sel(f_time=t_end)
        assert da_r.shape == da_y.shape, (
            f"{var} profile shape mismatch ref{da_r.shape} vs run{da_y.shape}"
        )
        x = np.asarray(da_r.values, dtype=np.float64).ravel()
        y = np.asarray(da_y.values, dtype=np.float64).ravel()
        ok = np.isclose(x, y, rtol=rtol, atol=atol, equal_nan=True)
        if not np.all(ok):
            bad = int(np.size(ok) - np.count_nonzero(ok))
            abs_err = np.abs(x - y)
            max_err = float(np.nanmax(abs_err))
            n_ref = ref.sizes["f_time"]
            n_run = run.sizes["f_time"]
            t_ref_last = np.asarray(ref["f_time"].isel(f_time=-1).values)
            t_run_last = np.asarray(run["f_time"].isel(f_time=-1).values)
            hint = (
                f"Diagnostics: {_env_digest_for_comparison()}; "
                f"f_time len ref={n_ref} run={n_run}; last ref={t_ref_last!r} run={t_run_last!r}. "
                "See module docstring: BLAS/threading/resample vs pinned reference NC."
            )
            raise AssertionError(
                f"{var}: {bad} / {ok.size} entries differ (last f_time slice); "
                f"max abs diff {max_err:g} (rtol={rtol}, atol={atol}). {hint}"
            )


@pytest.mark.slow
@pytest.mark.fortran
@pytest.mark.parametrize("var", ["T", "ALK", "NO3"])
def test_main_block_netcdf_last_timestep(main_block_paths, var):
    out_nc, _ = main_block_paths
    _array_last_time_close(var, _REF_NC, out_nc)
