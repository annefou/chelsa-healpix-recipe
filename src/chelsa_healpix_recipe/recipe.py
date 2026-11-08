"""The CHELSA recipe: source description, corrections, staging and healpix-convert settings."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from chelsa_healpix_recipe import conventions

#: The source, as cited. CHELSA v2.1 climatologies, 30 arc-seconds, EPSG:4326.
SOURCE = {
    "name": "CHELSA v2.1 climatologies",
    "doi": "10.16904/envidat.228",
    "license": "CC0-1.0",
    "version": "2.1",
    "reference_period": ["1981-01-01", "2010-12-31"],
    "specification": (
        "https://www.envidat.ch/dataset/2adb0c83-4653-4337-af28-f75c63ab7c74/resource/"
        "61fff0a7-6abb-45f7-a5f6-48f6e2c24851/download/chelsa_file_specification.pdf"
    ),
    "url_template": (
        "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/"
        "CHELSA_{name}_1981-2010_V.2.1.tif"
    ),
}

#: Variables with verified meaning, after the file's own scale and offset are applied.
#: The GeoTIFFs carry a scale and offset but no unit and no nodata value; units below
#: correct CHELSA's file specification where it is wrong (checked against CHELSA's own
#: monthly layers, see README).
VARIABLES = {
    "bio1": {"units": "degC", "standard_name": "air_temperature",
             "long_name": "Mean annual near-surface air temperature (BIO1)"},
    # specification says degC; the values are 100 x the population standard deviation
    "bio4": {"units": "0.01 degC",
             "long_name": "Temperature seasonality: 100 x standard deviation of monthly mean "
                          "temperatures (BIO4)"},
    "bio12": {"units": "kg m-2", "standard_name": "precipitation_amount",
              "long_name": "Annual precipitation amount (BIO12)"},
    # specification says kg m-2; the values are a coefficient of variation in percent
    "bio15": {"units": "percent",
              "long_name": "Precipitation seasonality: coefficient of variation of monthly "
                           "precipitation (BIO15)"},
}

_WGS84 = {"name": "wgs84"}


def settings(variable: str, depth: int = 13, pyramid_to: int = 8, chunk_depth: int = 7,
             convention: str = "as-published") -> dict:
    """healpix-convert settings for one variable, as plain data.

    CHELSA (about 1 km) is resampled bilinearly at ``depth`` (13 is about 800 m), then the
    pyramid builds the coarser depths down to ``pyramid_to`` with ``mean``: on equal-area
    cells, block means are area-weighted cell means.
    """
    meta = {**VARIABLES[variable], **conventions.attrs(variable, convention)}
    return {
        "metadata": {
            "stac": {"override_properties": {"license": SOURCE["license"]}},
        },
        "multiscale_settings": {
            "groups": ["chelsa"],
            "pyramid_groups": {
                "chelsa": {"derived_levels": list(range(depth - 1, pyramid_to - 1, -1)), "reduce": "mean"},
            },
        },
        "group_settings": {
            "chelsa/v21": {
                "healpix": {"refinement_level": depth, "indexing_scheme": "nested", "ellipsoid": _WGS84},
                "chunk": {
                    "method": "healpix_cell_dense",
                    "healpix": {"refinement_level": chunk_depth, "indexing_scheme": "nested", "ellipsoid": _WGS84},
                    "chunk_buffer_width": 3000.0,  # metres: a few CHELSA pixels
                },
                "resampler": {"name": "bilinear"},
                "output_dtype": "float32",
                "metadata": {"override_attrs": {variable: dict(meta)}},
            },
        },
    }


def _read_window(url: str, bbox: tuple[float, float, float, float], pad_deg: float):
    """CHELSA pixels in ``bbox`` + ``pad_deg``, scale and offset applied, rows south to north."""
    import rasterio
    from rasterio.windows import from_bounds

    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
    w, s, e, n = bbox[0] - pad_deg, bbox[1] - pad_deg, bbox[2] + pad_deg, bbox[3] + pad_deg
    with rasterio.open("/vsicurl/" + url) as src:
        win = from_bounds(w, s, e, n, src.transform).round_offsets().round_lengths()
        raw = src.read(1, window=win)
        scale, offset, nodata = src.scales[0], src.offsets[0], src.nodata
        t = src.window_transform(win)
    values = raw.astype("float64")
    if nodata is not None:  # CHELSA declares none (2026-11-04); kept for other sources
        values[raw == nodata] = np.nan
    values = values * scale + offset
    lon = t.c + t.a * (np.arange(values.shape[1]) + 0.5)
    lat = t.f + t.e * (np.arange(values.shape[0]) + 0.5)
    return values[::-1], lon, lat[::-1]  # ascending latitude


def stage(variable: str, bbox: tuple[float, float, float, float], path: str | Path,
          pad_deg: float = 0.5, convention: str = "as-published") -> Path:
    """Write the staged input healpix-convert reads: a lat/lon Zarr dataset with a STAC item."""
    import shapely
    import xarray as xr

    conventions.check(variable, convention)
    values, lon, lat = _read_window(SOURCE["url_template"].format(name=variable), bbox, pad_deg)
    bio12 = None
    if conventions.needs_bio12(variable, convention):
        bio12, _, _ = _read_window(SOURCE["url_template"].format(name="bio12"), bbox, pad_deg)
    values = conventions.apply(variable, values, convention, bio12=bio12)
    ds = xr.Dataset({variable: (("lat", "lon"), values.astype("float32"))}, coords={"lat": lat, "lon": lon})
    ds["lat"].attrs.update(standard_name="latitude", units="degrees_north")
    ds["lon"].attrs.update(standard_name="longitude", units="degrees_east")
    box = shapely.box(float(lon[0]), float(lat[0]), float(lon[-1]), float(lat[-1]))
    stac_item = {
        "type": "Feature",
        "stac_version": "1.1.0",
        # healpix-convert requires the field; the DEM GLO-30 notebook uses the projection extension
        "stac_extensions": ["https://stac-extensions.github.io/projection/v2.0.0/schema.json"],
        "id": f"chelsa-v21-{variable}-{conventions.target(convention)}",
        "geometry": shapely.geometry.mapping(box),
        "bbox": list(box.bounds),
        "properties": {
            "datetime": None,
            "start_datetime": "1981-01-01T00:00:00Z",
            "end_datetime": "2010-12-31T23:59:59Z",
            "proj:code": "EPSG:4326",
        },
        "links": [{"rel": "cite-as", "href": f"https://doi.org/{SOURCE['doi']}"}],
        "assets": {},
    }
    path = Path(path)
    xr.DataTree.from_dict({"/": xr.Dataset(attrs={"stac_discovery": stac_item}), "chelsa/v21": ds}).to_zarr(
        path, mode="w-", zarr_format=3)
    return path


def convert(variable: str, bbox: tuple[float, float, float, float], output: str | Path,
            staged: str | Path | None = None, convention: str = "as-published", **settings_kwargs):
    """Stage, then convert with healpix-convert's public API. Returns the output DataTree."""
    import healpix_convert as hc
    import shapely

    output = Path(output)
    staged = Path(staged) if staged else output.with_name(output.stem + "-staged.zarr")
    stage(variable, bbox, staged, convention=convention)
    cfg = hc.ConvertSettings.model_validate(settings(variable, convention=convention, **settings_kwargs))
    return hc.create_healpix_dataset([str(staged)], cfg, str(output), output_extent=shapely.box(*bbox))
