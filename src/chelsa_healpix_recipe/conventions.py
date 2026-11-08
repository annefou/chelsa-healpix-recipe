"""BIOCLIM conventions: the same variable names are computed differently by different sources.

- ``"population"``: population standard deviation (divide by n), CV of monthly precipitation
  (CHELSA v2.1).
- ``"dismo"``: sample standard deviation (divide by n-1), CV of monthly precipitation + 1 mm
  (R's dismo::biovars; WorldClim 2.1).

``convention="as-published"`` keeps the source's own values. Another convention is applied exactly:

- bio4 is rescaled by sqrt(11/12) or sqrt(12/11): both use the same 12 monthly values;
- bio15 also needs bio12 of the same source: adding 1 mm to every month leaves the standard
  deviation unchanged and raises the mean m = bio12 / 12 by 1, so
  CV_population = CV_dismo * (m + 1) / m * sqrt(11/12), and the inverse. Where m = 0 (no rain at
  all) the population CV is undefined (missing) and the dismo CV is exactly 0.

Checked against CV recomputed from the monthly layers: ratio 1.00000 (WorldClim to population) and
1.0000 +- 0.0009 (CHELSA to dismo; CHELSA stores values to 0.1). Variables defined the same way in
both conventions (bio1, bio12) are unchanged.
"""

from __future__ import annotations

import math

#: The convention this source uses.
SOURCE_CONVENTION = "population"  # CHELSA v2.1
CONVENTIONS = ("as-published", "population", "dismo")

_SD = {"population": "population standard deviation (n)", "dismo": "sample standard deviation (n-1)"}
_FACTOR = {("dismo", "population"): math.sqrt(11 / 12), ("population", "dismo"): math.sqrt(12 / 11)}
_BIO15 = {"population": "monthly precipitation, population standard deviation (n)",
          "dismo": "monthly precipitation plus 1 mm, sample standard deviation (n-1)"}


def target(convention: str) -> str:
    """The convention the output will be in."""
    if convention not in CONVENTIONS:
        raise ValueError(f"convention must be one of {CONVENTIONS}, got {convention!r}")
    return SOURCE_CONVENTION if convention == "as-published" else convention


def needs_bio12(variable: str, convention: str) -> bool:
    """Whether converting this variable needs bio12 of the same source and pixels."""
    return variable == "bio15" and target(convention) != SOURCE_CONVENTION


def apply(variable: str, values, convention: str, bio12=None):
    """The values in the requested convention (unchanged when it is the source's own).

    ``bio12`` (same source, same pixels, mm) is required to convert bio15.
    """
    import numpy as np

    to = target(convention)
    if to == SOURCE_CONVENTION:
        return values
    factor = _FACTOR[(SOURCE_CONVENTION, to)]
    if variable == "bio4":
        return values * factor
    if variable == "bio15":
        if bio12 is None:
            raise ValueError("converting bio15 needs bio12 of the same source and pixels")
        m = np.asarray(bio12, dtype="float64") / 12.0
        with np.errstate(divide="ignore", invalid="ignore"):
            if to == "population":  # from dismo: remove the +1 mm, then n-1 -> n
                return np.where(m > 0, values * (m + 1.0) / m * factor, np.nan)
            # from population: add the +1 mm, then n -> n-1. Where m = 0 every month is dry, so the
            # +1 mm months are all equal and the dismo CV is exactly 0 (the population CV was 0/0)
            return np.where(m > 0, values * m / (m + 1.0) * factor, 0.0)
    return values


def attrs(variable: str, convention: str) -> dict:
    """Attributes stating the convention of the output, and any conversion made."""
    to = target(convention)
    if variable == "bio4":
        out = {"long_name": f"Temperature seasonality: 100 x {_SD[to]} of monthly mean temperatures (BIO4)",
               "bioclim_convention": to}
        if to != SOURCE_CONVENTION:
            out["comment"] = (f"converted from the source's {_SD[SOURCE_CONVENTION]} to the {_SD[to]} "
                              f"by multiplying by {_FACTOR[(SOURCE_CONVENTION, to)]:.12f}")
        return out
    if variable == "bio15":
        out = {"long_name": f"Precipitation seasonality: coefficient of variation of {_BIO15[to]} (BIO15)",
               "bioclim_convention": to}
        if to != SOURCE_CONVENTION:
            out["comment"] = (f"converted from the source's convention ({_BIO15[SOURCE_CONVENTION]}) using "
                              "bio12 of the same source: m = bio12 / 12, CV_population = CV_dismo * (m + 1) / m "
                              "* sqrt(11/12) and its inverse; where m = 0, missing (population) or 0 (dismo)")
        return out
    return {}


def check(variable: str, convention: str) -> None:
    """Validate the convention name."""
    target(convention)
