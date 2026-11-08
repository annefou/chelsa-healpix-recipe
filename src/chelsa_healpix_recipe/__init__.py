"""CHELSA v2.1 bioclimatic variables on HEALPix (WGS84, NESTED), as a healpix-convert recipe.

A recipe is a staging step (CHELSA's GeoTIFFs read into a lat/lon Zarr dataset) plus a
healpix-convert settings dict, run with healpix-convert's public API. Nothing is added to
healpix-convert itself. The layout is a trial, discussed in
https://github.com/GRID4EARTH/healpix-convert/issues/18.
"""

from healpix_recipe_tools.bioclim import CONVENTIONS
from chelsa_healpix_recipe.recipe import SOURCE_CONVENTION, SOURCE, VARIABLES, convert, settings, stage

__version__ = "0.1.0.dev0"
__all__ = ["CONVENTIONS", "SOURCE_CONVENTION", "SOURCE", "VARIABLES", "convert", "settings", "stage"]
