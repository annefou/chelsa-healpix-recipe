"""Offline: the settings are valid healpix-convert settings and carry the corrections."""

import pytest

import healpix_convert as hc

from chelsa_healpix_recipe import VARIABLES, settings


@pytest.mark.parametrize("variable", sorted(VARIABLES))
def test_settings_validate_with_healpix_convert(variable):
    cfg = hc.ConvertSettings.model_validate(settings(variable))
    assert cfg is not None


@pytest.mark.parametrize("variable", sorted(VARIABLES))
def test_units_are_set_explicitly(variable):
    group = settings(variable)["group_settings"]["chelsa/v21"]
    assert group["metadata"]["override_attrs"][variable]["units"] == VARIABLES[variable]["units"]


def test_the_corrected_units_are_not_the_specification_ones():
    # CHELSA's file specification says degC for bio4 and kg m-2 for bio15
    assert VARIABLES["bio4"]["units"] == "0.01 degC"
    assert VARIABLES["bio15"]["units"] == "percent"


def test_the_pyramid_goes_from_the_converted_depth_down():
    s = settings("bio1", depth=13, pyramid_to=8)
    assert s["group_settings"]["chelsa/v21"]["healpix"]["refinement_level"] == 13
    assert s["multiscale_settings"]["pyramid_groups"]["chelsa"]["derived_levels"] == [12, 11, 10, 9, 8]
