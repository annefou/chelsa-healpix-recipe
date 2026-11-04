# Reference cell means

`bio1_iberia_cell_means.csv`: CHELSA v2.1 `bio1` (degC) averaged into HEALPix cells (WGS84,
NESTED) at depths 8 and 9, for the window 6-3W, 39-42N read with a 0.5 degree margin. Columns:
`depth`, `cell`, `mean_degC`, `coverage` (valid pixel area / cell area).

These values come from an independent implementation, not from this recipe: the area-weighted
binning of pixel centres in [healpix-connector](https://github.com/annefou/healpix-connector)
(`healpix_connector.binning.bin_to_cells`, commit `5a7e671`), written once on 2026-11-04. They
are stored here so that this repository does not depend on healpix-connector.

To regenerate them, in a healpix-connector environment:

```python
from healpix_connector.binning import bin_to_cells
from healpix_connector.sources import chelsa
from healpix_connector.sources.geotiff import read_window

win = read_window(chelsa.url(1), (-6.0, 39.0, -3.0, 42.0), pad_deg=0.5)
for depth in (8, 9):
    stats = bin_to_cells(win.values, win.lon, win.lat, depth, win.res_deg)
    # write depth, stats.cell_ids, stats.mean, stats.coverage
```
