"""Network: the converted cell means agree with healpix-connector's binning of the same pixels.

Run with ``pytest -m network``. Reads a 3 x 3 degree window of CHELSA bio1 over central Iberia.
Only coarse cells fully covered by data are compared: the pyramid's mean ignores missing
children, so partly covered cells look complete (see README).
"""

import numpy as np
import pytest
import xarray as xr

from chelsa_healpix_recipe import convert
from chelsa_healpix_recipe.recipe import SOURCE, _read_window

BBOX = (-6.0, 39.0, -3.0, 42.0)


@pytest.mark.network
@pytest.mark.parametrize("depth, mean_tol, max_tol", [(8, 0.005, 0.03), (9, 0.01, 0.1)])
def test_pyramid_means_match_healpix_connector(tmp_path, depth, mean_tol, max_tol):
    from healpix_connector.binning import bin_to_cells

    convert("bio1", BBOX, tmp_path / "out.zarr")
    dt = xr.open_datatree(tmp_path / "out.zarr", engine="zarr")
    out = dt[f"chelsa/{depth}"].to_dataset()
    got = dict(zip(out["cell_ids"].values.astype(np.uint64).tolist(), out["bio1"].values.astype(float)))

    values, lon, lat = _read_window(SOURCE["url_template"].format(name="bio1"), BBOX, 0.5)
    ref = bin_to_cells(values, lon, lat, depth, 1 / 120)
    pairs = [(got[int(c)], m) for c, m, cov in zip(ref.cell_ids, ref.mean, ref.coverage)
             if cov >= 0.999 and int(c) in got and np.isfinite(got[int(c)])]
    assert len(pairs) > 50
    diff = np.array([a - b for a, b in pairs])
    assert np.abs(diff).mean() < mean_tol
    assert np.abs(diff).max() < max_tol
    assert abs(diff.mean()) < mean_tol  # no bias


@pytest.mark.network
def test_units_reach_the_output(tmp_path):
    convert("bio4", (-5.0, 40.0, -4.0, 41.0), tmp_path / "bio4.zarr", pyramid_to=11)
    dt = xr.open_datatree(tmp_path / "bio4.zarr", engine="zarr")
    assert dt["chelsa/13"].to_dataset()["bio4"].attrs["units"] == "0.01 degC"
