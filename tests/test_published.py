"""Counts published by Hart & Forcade, Phys. Rev. B 77, 224115 (2008)
(all compositions, every species present, label exchange applied)."""
import pytest

import necklace_enum as te
from necklace_enum.lattice import unique_superlattices

# Table V: binary fcc, n = 2..23
FCC2 = [2, 3, 12, 14, 50, 52, 229, 252, 685, 682, 3875, 2624, 9628, 16584, 49764, 42135,
        212612, 174104, 867893, 1120708, 2628180, 3042732]
# Table VII: ternary and quaternary fcc
FCC3 = {3: 3, 4: 13, 5: 23, 6: 130, 7: 197, 8: 1267, 9: 2322, 10: 9332}
FCC4 = {4: 7, 5: 9, 6: 110, 7: 211, 8: 2110, 9: 5471, 10: 32362}
# Table VI: simple cubic
SC = {2: 3, 3: 3, 4: 15}
# Table IV: fcc superlattices (OEIS A159842)
FCC_SUPER = [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]


def test_fcc_superlattices():
    rots = te.fcc().rots
    assert [len(unique_superlattices(n, rots)) for n in range(1, 11)] == FCC_SUPER


def test_fcc_binary():
    got = te.count_structures(te.fcc(2), range(2, 18), label_exchange=True)
    assert [got[n] for n in range(2, 18)] == FCC2[:16]


@pytest.mark.slow
def test_fcc_binary_large():
    got = te.count_structures(te.fcc(2), range(18, 24), label_exchange=True)
    assert [got[n] for n in range(18, 24)] == FCC2[16:]


def test_fcc_ternary():
    got = te.count_structures(te.fcc(3), list(FCC3), label_exchange=True)
    assert got == FCC3


def test_fcc_quaternary():
    got = te.count_structures(te.fcc(4), list(FCC4), label_exchange=True)
    assert got == FCC4


def test_simple_cubic():
    assert te.count_structures(te.sc(2), list(SC), label_exchange=True) == SC
