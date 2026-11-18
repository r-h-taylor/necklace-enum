"""Exported crystals are correct, primitive in the expected cell and pairwise distinct."""
import itertools

import numpy as np
import pytest

import necklace_enum as te

ase = pytest.importorskip("ase")
spglib = pytest.importorskip("spglib")
from ase.utils.structure_comparator import SymmetryEquivalenceCheck  # noqa: E402


def _primitive_atoms(atoms):
    cell = (atoms.cell[:], atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    prim = spglib.find_primitive(cell, symprec=1e-5)
    return len(prim[2])


@pytest.mark.parametrize("parent,sizes", [(te.fcc(["Cu", "Au"]), range(1, 8)),
                                          (te.hcp(["Mg", "Cd"]), range(1, 4))])
def test_export_distinct_and_primitive(parent, sizes):
    structs = list(te.enumerate_structures(parent, sizes))
    atoms = [s.to_ase() for s in structs]
    for s, a in zip(structs, atoms):
        assert len(a) == s.size * parent.d
        assert _primitive_atoms(a) == len(a), "structure is periodic in a smaller cell"
    comp = SymmetryEquivalenceCheck()
    # candidates must share formula and space group; ASE's comparator alone can report
    # false matches between structures of different symmetry
    by_formula = {}
    for a in atoms:
        sg = spglib.get_symmetry_dataset((a.cell[:], a.get_scaled_positions(), a.get_atomic_numbers()),
                                         symprec=1e-5).number
        by_formula.setdefault((a.get_chemical_formula(), sg), []).append(a)
    for group in by_formula.values():
        for x, y in itertools.combinations(group, 2):
            assert not comp.compare(x, y), "two exported structures are equivalent"


def test_composition_dict():
    P = te.fcc(["Cu", "Au"])
    structs = list(te.enumerate_structures(P, 8, composition={"Cu": 6, "Au": 2}))
    assert len(structs) == te.count_structures(P, 8, composition=(6, 2))[8]
    assert all(s.to_ase().get_chemical_formula() == "Au2Cu6" for s in structs)
