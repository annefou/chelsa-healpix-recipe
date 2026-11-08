# chelsa-healpix-recipe

[CHELSA v2.1](https://doi.org/10.16904/envidat.228) bioclimatic variables on HEALPix
(WGS84 ellipsoid, NESTED cell ids), converted with
[healpix-convert](https://github.com/GRID4EARTH/healpix-convert) as a third-party recipe.

> **Status: trial.** This is a first try at "one repository per dataset recipe". The layout is
> being discussed in [healpix-convert#18](https://github.com/GRID4EARTH/healpix-convert/issues/18)
> and may change.

## What a recipe is here

Nothing is added to healpix-convert. A recipe is:

1. **a staging step**: CHELSA's GeoTIFFs are read with rasterio (the file's own scale and offset
   applied) and written as a lat/lon Zarr dataset with a STAC item, which is what
   healpix-convert reads, following the DEM GLO-30 example of healpix-convert PR #76;
2. **a settings dict**, plain data validated by healpix-convert: bilinear resampling to depth 13
   (about 800 m, close to CHELSA's 1 km), then a pyramid of `mean` down to depth 8. On
   equal-area cells, block means are area-weighted cell means;
3. **corrections**: units set explicitly with `override_attrs`, because CHELSA's GeoTIFFs carry
   no unit and its file specification is wrong for two variables (below).

```python
from chelsa_healpix_recipe import convert

dt = convert("bio1", bbox=(-6.0, 39.0, -3.0, 42.0), output="bio1-iberia.zarr")
dt["chelsa/8"]      # depth 8 (about 25 km) cell means
```

## Variables

| Variable | Unit | Note |
|---|---|---|
| `bio1` | degC | mean annual near-surface air temperature |
| `bio4` | 0.01 degC | 100 x the population standard deviation of monthly means; the specification says degC |
| `bio12` | kg m-2 | annual precipitation |
| `bio15` | percent | coefficient of variation of monthly precipitation; the specification says kg m-2 |

The `bio4` and `bio15` corrections were checked against CHELSA's own monthly layers (a 1 degree
window at 39N 6W): ratio 100.00 to the population standard deviation, and 1.000 to the
coefficient of variation. See healpix-connector's
[`examples/verify_chelsa_scaling.py`](https://github.com/annefou/healpix-connector/blob/main/examples/verify_chelsa_scaling.py).

## Conventions: comparing with WorldClim

CHELSA and WorldClim give the same BIOCLIM names to differently computed values: CHELSA uses the
population standard deviation (n) and the CV of monthly precipitation; WorldClim follows R's
`dismo::biovars`, with the sample standard deviation (n-1) and the CV of monthly precipitation
plus 1 mm. `bio4` differs by sqrt(12/11), about 4.4 %, everywhere.

`convention=` converts exactly, so values from the two sources can be compared:

```python
convert("bio4",  bbox, "bio4.zarr",  convention="dismo")    # as WorldClim computes it
convert("bio15", bbox, "bio15.zarr", convention="dismo")    # also reads bio12 (see below)
```

| `convention` | Meaning |
|---|---|
| `"as-published"` (default) | the source's own values (population convention) |
| `"population"` | population SD, CV of monthly precipitation (CHELSA's) |
| `"dismo"` | sample SD, CV of monthly precipitation + 1 mm (WorldClim's) |

- `bio4`: rescaled by sqrt(11/12) or sqrt(12/11), exact (both use the same 12 monthly values).
- `bio15`: adding 1 mm to every month leaves the standard deviation unchanged and raises the mean
  m = bio12 / 12 by 1, so the conversion uses `bio12` of the same source:
  CV_population = CV_dismo x (m + 1) / m x sqrt(11/12), and its inverse. Where it never rains
  (m = 0) the population CV is missing and the dismo CV is 0. Checked against the CV recomputed
  from the monthly layers: ratio 1.00000 (WorldClim to population) and 1.0000 +- 0.0009 (CHELSA to
  dismo, whose values are stored to 0.1).
- `bio1`, `bio12`: defined alike, unchanged.

The output states its convention (`bioclim_convention`) and, when converted, how (`comment`).

## Agreement with an independent reference

The network tests convert `bio1` over central Iberia (6-3W, 39-42N) and compare the pyramid
means with stored reference values computed independently (an area-weighted binning of the same
pixels, see [tests/reference](tests/reference/README.md)), on cells fully covered by data.
This repository does not depend on the code that produced them. Measured 2026-11-04
(healpix-convert `dfa4ed3`):

| Depth | Cells | Mean difference | Mean abs. difference | Max abs. difference |
|---|---|---|---|---|
| 8 (~25 km) | 68 | 0.000 degC | 0.0015 degC | 0.010 degC |
| 9 (~12 km) | 366 | 0.000 degC | 0.0038 degC | 0.054 degC |

The conversion of the window took about 5 s (CPU only).

## Open points

- **Coverage is not reported.** The pyramid's `mean` ignores missing children, so a coarse cell
  only partly covered by data (a window edge, a coastline, a gap) gets a value that describes
  only part of its area, with nothing saying so. A coverage fraction, and the spread within each
  cell, would need reductions beyond `mean` in healpix-convert (or
  [healpix-resample#63](https://github.com/GRID4EARTH/healpix-resample/issues/63)).
- **CHELSA declares no nodata value.** How areas without data (ocean, gaps) are encoded is still
  to be checked; if they hold a number, the conversion averages it in. The trial window is
  inland, so it did not show.
- **The staged STAC item needs `stac_extensions`**; healpix-convert fails without it.

## Running

```bash
pixi run -e test test            # offline: the settings are valid healpix-convert settings
pixi run -e test test-network    # reads CHELSA and compares with the stored reference
```

## Licence

Code: Apache License 2.0. CHELSA data: CC0-1.0; cite Karger et al., CHELSA v2.1,
[doi:10.16904/envidat.228](https://doi.org/10.16904/envidat.228), and Karger et al. (2017),
*Scientific Data* 4, 170122, [doi:10.1038/sdata.2017.122](https://doi.org/10.1038/sdata.2017.122).
