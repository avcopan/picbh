"""Graph tests."""

import enum

import pytest

from picbh import graph
from picbh.graph.base import FieldEnumMismatchError, FieldEnumModel


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
        ("C", []),
        ("CC", [1]),
        ("C=C", [2]),
        ("C#C", [3]),
        ("CC=O", [1, 2]),
    ],
)
def test__bond_orders(smi: str, ref: list[int]) -> None:
    """Test graph bond orders."""
    gra = graph.from_smiles(smi)
    assert graph.bond_orders(gra) == ref


def test__bond_orders__keys() -> None:
    """Bond orders are returned for the requested bonds, in the requested order."""
    gra = graph.from_smiles("C=CC#N")  # bonds: (0, 1) = 2, (1, 2) = 1, (2, 3) = 3
    assert graph.bond_orders(gra, [(2, 3), (0, 1)]) == [3, 2]


def test__bond_orders__kekulized_aromatic() -> None:
    """Aromatic bonds are Kekulized to alternating single and double orders."""
    benzene = graph.from_smiles("c1ccccc1")
    assert sorted(graph.bond_orders(benzene)) == [1, 1, 1, 2, 2, 2]


def test__subgraph() -> None:
    """The induced subgraph keeps only the given nodes and their mutual bonds."""
    propane = graph.from_smiles("CCC")
    ethyl = graph.subgraph(propane, [0, 1])
    assert graph.is_isomorphic(ethyl, graph.from_smiles("[CH3][CH2]"))
    assert graph.node_keys(propane) == [0, 1, 2]  # copy left the input untouched


def test__subgraph__in_place() -> None:
    """An in-place subgraph mutates and returns the original graph."""
    propane = graph.from_smiles("CCC")
    ethyl = graph.subgraph(propane, [0, 1], in_place=True)
    assert ethyl is propane
    assert graph.node_keys(propane) == [0, 1]


def test__neighborhood() -> None:
    """The neighborhood collects nodes within a radius of any center node."""
    chain = graph.from_smiles("CCCCC")  # 0 - 1 - 2 - 3 - 4
    assert graph.neighborhood(chain, 2, radius=1) == frozenset({1, 2, 3})
    assert graph.neighborhood(chain, (0, 4), radius=1) == frozenset({0, 1, 3, 4})


def test__capped_subgraph() -> None:
    """A severed bond is capped with implicit hydrogens equal to its bond order."""
    propene = graph.from_smiles("C=CC")  # 0 = 1 - 2
    frag = graph.capped_subgraph(propene, [1, 2])
    assert graph.inchi(frag) == graph.inchi(graph.from_smiles("CC"))


def test__field_enum_matches_fields() -> None:
    """A nested Field enum in parity with the model fields is accepted."""

    class Good(FieldEnumModel):
        class Field(enum.StrEnum):
            a = "a"
            b = "b"

        a: int
        b: str

    assert Good.Field.a == "a"
    assert set(Good.Field) == set(Good.model_fields)


def test__field_enum_missing_member() -> None:
    """A field with no Field member fails at class-creation time."""
    with pytest.raises(FieldEnumMismatchError, match="missing from Field"):

        class MissingMember(FieldEnumModel):
            class Field(enum.StrEnum):
                a = "a"

            a: int
            b: str


def test__field_enum_extra_member() -> None:
    """A Field member that is not a field fails at class-creation time."""
    with pytest.raises(FieldEnumMismatchError, match="not a model field"):

        class ExtraMember(FieldEnumModel):
            class Field(enum.StrEnum):
                a = "a"
                c = "c"

            a: int


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
