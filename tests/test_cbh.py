"""CBH tests."""

from collections import Counter, defaultdict

import pytest

from picbh import cbh, graph
from picbh.graph import Atom

# Methionine, the worked example of Scheme 1 in the CBH automation paper
# (Sengupta, Raghavachari et al., J. Chem. Theory Comput. 2023, 19, 3763).
METHIONINE = "CSCCC(N)C(=O)O"

# Each rung as drawn in Scheme 1, as {SMILES: coefficient} with S = Sum_F nu(F) F
# (positive = product side, negative = reactant side).
SCHEME_1: dict[int, dict[str, int]] = {
    0: {"C": 5, "N": 1, "O": 2, "S": 1, "[H][H]": -9},
    1: {"CC": 3, "CN": 1, "CS": 2, "C=O": 1, "CO": 1, "C": -6, "S": -1},
    2: {
        "CSC": 1,
        "CCS": 1,
        "CCC": 1,
        "CC(C)N": 1,
        "CC(=O)O": 1,
        "CC": -3,
        "CS": -1,
    },
    3: {
        "CCSC": 1,
        "CCCS": 1,
        "CCC(C)N": 1,
        "CC(N)C(=O)O": 1,
        "CCC": -1,
        "CCS": -1,
        "CC(C)N": -1,
    },
    4: {
        "CCCSC": 1,
        "CC(N)CCS": 1,
        "CCC(N)C(=O)O": 1,
        "CCCS": -1,
        "CCC(C)N": -1,
    },
}


def _by_inchi(spec: dict[str, int]) -> dict[str, int]:
    """Convert a {SMILES: coefficient} spec to {InChI: coefficient}, dropping zeros."""
    coeffs: dict[str, int] = defaultdict(int)
    for smi, coeff in spec.items():
        coeffs[graph.inchi(graph.from_smiles(smi))] += coeff
    return {frag: coeff for frag, coeff in coeffs.items() if coeff}


def _element_counts(gra: graph.MolGraph) -> Counter[str]:
    """Count every atom in a molecular graph, including implicit hydrogens."""
    counts = Counter(graph.symbols(gra))
    counts["H"] += sum(
        gra.nodes[key][Atom.Field.implicit_hydrogens] for key in graph.node_keys(gra)
    )
    return counts


@pytest.mark.parametrize("rung", [0, 1, 2, 3, 4])
def test__expansion__methionine(rung: int) -> None:
    """CBH(0-4) for methionine reproduces Scheme 1 of the CBH automation paper."""
    gra = graph.from_smiles(METHIONINE)
    assert cbh.expansion(gra, rung=rung) == _by_inchi(SCHEME_1[rung])


@pytest.mark.parametrize("rung", [0, 1, 2, 3, 4])
def test__expansion__stoichiometry_balanced(rung: int) -> None:
    """Every CBH rung conserves each element, including hydrogen."""
    gra = graph.from_smiles(METHIONINE)
    total: Counter[str] = Counter()
    for chi, coeff in cbh.expansion(gra, rung=rung).items():
        for element, count in _element_counts(graph.from_inchi(chi)).items():
            total[element] += coeff * count
    assert +total == _element_counts(gra)
