"""``tumult.linelist.search_lines``: the line-list panel's free-text search."""
import pandas as pd
import pytest

from tumult.linelist import search_lines


@pytest.fixture
def lines():
    return pd.DataFrame({
        "ion": ["O I", "O I", "O I", "O II", "OH", "Os I", "Fe I", "Fe II", "CH"],
        "lam0": [7771.959, 7774.177, 7775.398, 7772.5, 3100.1, 4420.5, 5000.47, 5018.44,
                 43000.0],
    })


def found(df):
    return list(zip(df.ion, df.lam0))


def test_species_words_narrow_as_they_are_typed(lines):
    assert set(search_lines(lines, "O").ion) == {"O I", "O II", "OH", "Os I"}
    assert set(search_lines(lines, "O I").ion) == {"O I", "O II"}
    assert set(search_lines(lines, "O I ").ion) == {"O I"}
    assert set(search_lines(lines, "fe ii").ion) == {"Fe II"}
    assert set(search_lines(lines, "Fe").ion) == {"Fe I", "Fe II"}


def test_digits_are_the_leading_digits_of_the_wavelength(lines):
    assert found(search_lines(lines, "O I 77")) == [
        ("O I", 7771.959), ("O I", 7774.177), ("O I", 7775.398)]
    assert found(search_lines(lines, "O I 777")) == found(search_lines(lines, "O I 77"))
    assert found(search_lines(lines, "O I 7774")) == [("O I", 7774.177)]
    # a number after the species completes it: O II is not O I any more
    assert "O II" not in set(search_lines(lines, "O I 77").ion)


def test_wavelength_reads_rounded_or_truncated(lines):
    assert found(search_lines(lines, "7772")) == [("O I", 7771.959), ("O II", 7772.5)]
    assert found(search_lines(lines, "7771")) == [("O I", 7771.959)]
    assert found(search_lines(lines, "7771.9")) == [("O I", 7771.959)]
    assert found(search_lines(lines, "7771.96")) == [("O I", 7771.959)]
    assert found(search_lines(lines, "7771.95")) == [("O I", 7771.959)]  # truncated
    assert found(search_lines(lines, "7771.94")) == []


def test_decimal_point_pins_the_magnitude_and_integers_do_not(lines):
    assert found(search_lines(lines, "43")) == [("CH", 43000.0)]
    assert found(search_lines(lines, "43.")) == []
    assert found(search_lines(lines, "50")) == [("Fe I", 5000.47), ("Fe II", 5018.44)]


def test_several_wavelengths_match_any(lines):
    assert found(search_lines(lines, "O I 7772 7775")) == [
        ("O I", 7771.959), ("O I", 7775.398)]


def test_empty_query_matches_nothing(lines):
    assert len(search_lines(lines, "")) == 0
    assert len(search_lines(lines, "   ")) == 0
