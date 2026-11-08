"""The CHELSA recipe: what is specific to CHELSA v2.1. The rest comes from healpix-recipe-tools."""

from __future__ import annotations

import dataclasses
from pathlib import Path

from healpix_recipe_tools import (DEFAULT_LAND_MASK, apply_land_mask, bioclim, gdal_path, mask_attrs,
                                  raster_settings, read_window, run, stac_item, write_staged)

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

#: CHELSA's BIOCLIM convention: population standard deviation, CV of monthly precipitation.
SOURCE_CONVENTION = bioclim.SOURCES["CHELSA v2.1"]

GROUP = "chelsa/v21"

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


def settings(variable: str, depth: int = 13, pyramid_to: int = 8, chunk_depth: int = 7,
             convention: str = "as-published", land_mask=DEFAULT_LAND_MASK) -> dict:
    """healpix-convert settings for one variable, as plain data (see healpix_recipe_tools.raster_settings)."""
    attrs = {**VARIABLES[variable], **bioclim.attrs(variable, convention, SOURCE_CONVENTION), **mask_attrs(land_mask)}
    return raster_settings(GROUP, variable, attrs, depth=depth, pyramid_to=pyramid_to,
                           chunk_depth=chunk_depth, license=SOURCE["license"])


def _read(variable: str, bbox, pad_deg: float):
    return read_window(gdal_path(SOURCE["url_template"].format(name=variable)), bbox, pad_deg)


def stage(variable: str, bbox: tuple[float, float, float, float], path: str | Path,
          pad_deg: float = 0.5, convention: str = "as-published", land_mask=DEFAULT_LAND_MASK) -> Path:
    """Write the staged input healpix-convert reads: a lat/lon Zarr dataset with a STAC item.

    ``land_mask``: see healpix_recipe_tools.masks (default Natural Earth 10 m land; None keeps the sea).
    """
    bioclim.check(variable, convention, SOURCE_CONVENTION)
    window = _read(variable, bbox, pad_deg)
    bio12 = (_read("bio12", bbox, pad_deg).values
             if bioclim.needs_bio12(variable, convention, SOURCE_CONVENTION) else None)
    window = dataclasses.replace(
        window, values=bioclim.apply(variable, window.values, convention, SOURCE_CONVENTION, bio12=bio12))
    window = apply_land_mask(window, land_mask)  # the user's land mask (None keeps the sea)
    item = stac_item(
        f"chelsa-v21-{variable}-{bioclim.target(convention, SOURCE_CONVENTION)}", window,
        tuple(SOURCE["reference_period"]), license=SOURCE["license"],
        links=[{"rel": "cite-as", "href": f"https://doi.org/{SOURCE['doi']}"}])
    return write_staged(path, GROUP, variable, window, item)


def convert(variable: str, bbox: tuple[float, float, float, float], output: str | Path,
            staged: str | Path | None = None, convention: str = "as-published", land_mask=DEFAULT_LAND_MASK,
            **settings_kwargs):
    """Stage, then convert with healpix-convert's public API. Returns the output DataTree."""
    output = Path(output)
    staged = Path(staged) if staged else output.with_name(output.stem + "-staged.zarr")
    stage(variable, bbox, staged, convention=convention, land_mask=land_mask)
    return run(staged, settings(variable, convention=convention, land_mask=land_mask, **settings_kwargs),
               output, bbox)
