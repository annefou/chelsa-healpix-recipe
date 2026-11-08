"""Offline: the convention conversions are exact, checked on synthetic monthly data."""

import numpy as np
import pytest

import healpix_convert as hc

from chelsa_healpix_recipe import conventions, settings

rng = np.random.default_rng(0)
# 1000 synthetic places x 12 months, including dry months and a fully dry place
TAS = rng.normal(12, 6, (1000, 12))
PREC = rng.gamma(0.8, 40, (1000, 12))
PREC[:50, 5:8] = 0.0
PREC[0] = 0.0


def published(variable, convention):
    """What a source following ``convention`` publishes, computed from the monthly values."""
    ddof = 0 if convention == "population" else 1
    if variable == "bio4":
        return 100 * TAS.std(1, ddof=ddof)
    p = PREC if convention == "population" else PREC + 1
    with np.errstate(divide="ignore", invalid="ignore"):
        return 100 * p.std(1, ddof=ddof) / p.mean(1)


BIO12 = PREC.sum(1)
OTHER = "dismo" if conventions.SOURCE_CONVENTION == "population" else "population"


@pytest.mark.parametrize("variable", ["bio4", "bio15"])
def test_converting_to_the_other_convention_gives_its_published_values(variable):
    got = conventions.apply(variable, published(variable, conventions.SOURCE_CONVENTION), OTHER, bio12=BIO12)
    want = published(variable, OTHER)
    ok = np.isfinite(want)
    np.testing.assert_allclose(got[ok], want[ok], rtol=1e-12, atol=1e-12)
    assert np.isnan(got[~ok]).all()  # an undefined CV stays missing, never a number
    if variable == "bio15":  # the fully dry place: 0/0 in the population convention, exactly 0 in dismo
        assert (np.isnan(got[0]) if OTHER == "population" else got[0] == 0.0)


def test_as_published_changes_nothing():
    v = published("bio4", conventions.SOURCE_CONVENTION)
    assert conventions.apply("bio4", v, "as-published") is v
    assert conventions.apply("bio4", v, conventions.SOURCE_CONVENTION) is v


def test_variables_defined_alike_are_unchanged():
    v = TAS.mean(1)
    assert conventions.apply("bio1", v, OTHER) is v


def test_bio15_needs_bio12():
    with pytest.raises(ValueError, match="bio12"):
        conventions.apply("bio15", published("bio15", conventions.SOURCE_CONVENTION), OTHER)
    assert conventions.needs_bio12("bio15", OTHER) and not conventions.needs_bio12("bio4", OTHER)


def test_unknown_convention_is_refused():
    with pytest.raises(ValueError, match="convention"):
        conventions.target("anuclim")


@pytest.mark.parametrize("variable", ["bio4", "bio15"])
@pytest.mark.parametrize("convention", ["as-published", "population", "dismo"])
def test_the_output_says_its_convention(variable, convention):
    a = conventions.attrs(variable, convention)
    assert a["bioclim_convention"] == conventions.target(convention)
    assert ("comment" in a) == (conventions.target(convention) != conventions.SOURCE_CONVENTION)
    meta = settings(variable, convention=convention)
    assert hc.ConvertSettings.model_validate(meta) is not None
