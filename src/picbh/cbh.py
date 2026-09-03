"""Connectivity-based hierarchy (CBH).

Expands a target species into a stoichiometrically balanced set of smaller reference
species whose local connectivity matches the target's, via inclusion-exclusion over the
nerve of the primary fragments (see ``CBH Notes (Updated).typ``). Rung ``n`` preserves
every atom-centered (even ``n``) or bond-centered (odd ``n``) environment of radius
``n // 2``.
"""

import itertools
from collections import defaultdict

from .graph.base import neighborhood
from .graph.mol import MolGraph, bond_orders, capped_subgraph, inchi

# Standard InChI for dihydrogen, the CBH-0 reactant that hydrogenates each bond.
_H2_INCHI = "InChI=1S/H2/h1H"


def expansion(gra: MolGraph, *, rung: int) -> dict[str, int]:
    """Expand a species into its CBH(``rung``) reference species.

    Args:
        gra: A molecular graph of the target species.
        rung: The CBH rung.

    Returns:
        A mapping from canonical InChI to signed stoichiometric coefficient, such that
        ``S = Sum_F nu(F) * F``: positive coefficients are products, negative
        coefficients are reactants, and the target species itself is not included.
        Terms that cancel to zero are omitted.
    """
    fragments = primary_fragments(gra, rung=rung)
    coeffs: dict[str, int] = defaultdict(int)
    for simplex in _nerve(gra, fragments):
        nodes = frozenset.intersection(*(fragments[i] for i in simplex))
        coeffs[inchi(capped_subgraph(gra, nodes))] += (-1) ** (len(simplex) + 1)
    if rung == 0:
        coeffs[_H2_INCHI] -= sum(bond_orders(gra))
    return {frag: coeff for frag, coeff in coeffs.items() if coeff}


def primary_fragments(gra: MolGraph, *, rung: int) -> list[frozenset[int]]:
    """Get the node sets of the primary fragments for a CBH rung.

    A primary fragment whose atoms are contained in another one contributes only terms
    that cancel in `expansion` (its singleton is annulled by pairing every simplex that
    contains it with that simplex plus the containing fragment), so only the maximal
    fragments are kept. In a chain these are the centers at least ``rung // 2`` bonds
    from a terminus.

    Args:
        gra: A molecular graph.
        rung: The CBH rung.

    Returns:
        One node set per maximal primary center (atoms for even rungs, bonds for odd
        rungs), sorted so that a set's position is its fragment id.
    """
    radius = rung // 2
    centers = gra.edges if rung % 2 else ((key,) for key in gra.nodes)
    fragments = {neighborhood(gra, center, radius=radius) for center in centers}
    maximal = [
        frag for frag in fragments if not any(frag < other for other in fragments)
    ]
    return sorted(maximal, key=sorted)


def _nerve(gra: MolGraph, fragments: list[frozenset[int]]) -> set[frozenset[int]]:
    """Get the nerve of the primary fragments.

    The nerve ``N(X) = Union_v P+(X_v)`` is every set of fragment ids whose primary
    fragments share at least one atom; such a set maps to the overlap fragment
    ``H_Y = Intersection_{i in Y} H_i``. Since ``P+`` is monotonic, only the incidence
    sets ``X_v`` that are not contained in another one need to be expanded.

    Args:
        gra: A molecular graph.
        fragments: The primary fragment node sets, indexed by fragment id.

    Returns:
        The nerve, as sets of fragment ids.
    """
    incidences = {
        frozenset(i for i, frag in enumerate(fragments) if node in frag) for node in gra
    }
    maximal = [
        inc for inc in incidences if not any(inc < other for other in incidences)
    ]
    simplices: set[frozenset[int]] = set()
    for incidence in maximal:
        members = sorted(incidence)
        for size in range(1, len(members) + 1):
            simplices.update(map(frozenset, itertools.combinations(members, size)))
    return simplices
