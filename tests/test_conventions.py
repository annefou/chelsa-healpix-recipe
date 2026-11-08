"""Offline: this recipe declares its BIOCLIM convention and uses the shared, tested conversions.

The conversions themselves are tested in healpix-recipe-tools (tests/test_bioclim.py).
"""

import healpix_convert as hc
import pytest

from chelsa_healpix_recipe import CONVENTIONS, SOURCE_CONVENTION, settings


def test_the_source_convention_is_known():
    assert SOURCE_CONVENTION in ("population", "dismo")


@pytest.mark.parametrize("variable", ["bio4", "bio15"])
@pytest.mark.parametrize("convention", CONVENTIONS)
def test_every_convention_gives_valid_settings_that_state_it(variable, convention):
    s = settings(variable, convention=convention)
    group = next(iter(s["group_settings"].values()))
    attrs = group["metadata"]["override_attrs"][variable]
    expected = SOURCE_CONVENTION if convention == "as-published" else convention
    assert attrs["bioclim_convention"] == expected
    assert hc.ConvertSettings.model_validate(s) is not None
