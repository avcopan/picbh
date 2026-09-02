"""Graph tests."""

import pytest

from picbh import graph


def test__smiles() -> None:
    """Test graph smiles."""
    water_smiles = "O"
    water_inchi = "InChI=1S/H2O/h1H2"
    water = graph.from_smiles(water_smiles)
    assert graph.inchi(water) == water_inchi


def test__inchi() -> None:
    """Test graph inchi."""
    water_inchi = "InChI=1S/H2O/h1H2"
    water = graph.from_inchi(water_inchi)
    assert graph.inchi(water) == water_inchi


def test__remove_edges() -> None:
    """Test graph remove edges."""
    methanol = graph.from_smiles("CO")
    ch3_oh_ref = graph.from_smiles("[CH3].[OH]")
    ch3_oh = graph.remove_edges(methanol, [(0, 1)])
    assert graph.is_isomorphic(ch3_oh, ch3_oh_ref)


def test__symbols() -> None:
    """Test graph symbols."""
    water = graph.from_smiles("O")
    assert graph.symbols(water) == ["O"]


@pytest.mark.parametrize(
    ("smi", "ref"),
    [
        ("[H]", [1]),
        ("[He]", [0]),
        ("O", [0]),
        ("[OH]", [1]),
        ("[CH2]", [2]),
        ("C=C", [0, 0]),
        ("C#C", [0, 0]),
        ("O[O]", [0, 1]),
    ],
)
def test__unpaired_electrons(smi: str, ref: list[int]) -> None:
    """Test graph unpaired electrons."""
    gra = graph.from_smiles(smi)
    assert graph.unpaired_electrons(gra) == ref


@pytest.mark.parametrize(
    "smi",
    [
        "O",
        "[H]",
        "[H][H]",
        "[CH2]",
        "C=C",
        "C#C",
        "O[O]",
        "[OH]",
        "c1ccccc1",
        "c1ccc2ccccc2c1",
        "c1cc[nH]c1",
        "CC(=O)O",
        "FC(F)(F)F",
        "C1CC1",
        "N",
        "CN",
    ],
)
def test__smiles_roundtrip(smi: str) -> None:
    """Test SMILES <-> graph round trip.

    An unmodified graph converts back to the exact SMILES it was read from.
    """
    gra = graph.from_smiles(smi)
    assert graph.smiles(gra) == smi
