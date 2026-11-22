"""The Burnside counter against the generator and against published counts.

The two routes share only the superlattice list and the symmetry operations: the
generator never forms a group-element sum, and the Burnside counter never builds
a word, so agreement checks both.
"""
import itertools

import pytest

import necklace_enum as te
from necklace_enum import count_structures_burnside as burnside
from test_published import FCC2, FCC3, FCC4, SC

# binary fcc up to label exchange, n = 24..32, from the generator (scripts/reproduce_full_output.txt)
FCC2_LARGE = {24: 19335101, 25: 15592901, 26: 45995141, 27: 82823587, 28: 255215928,
              29: 224864930, 30: 1311471651, 31: 941857818, 32: 4592213143}


def _agree(P, sizes, **kw):
    assert burnside(P, sizes, **kw) == te.count_structures(P, sizes, **kw)


def test_published_hart_forcade():
    got = burnside(te.fcc(2), range(2, 24), label_exchange=True)
    assert [got[n] for n in range(2, 24)] == FCC2
    assert burnside(te.fcc(3), list(FCC3), label_exchange=True) == FCC3
    assert burnside(te.fcc(4), list(FCC4), label_exchange=True) == FCC4
    assert burnside(te.sc(2), list(SC), label_exchange=True) == SC


def test_binary_fcc_beyond_published():
    assert burnside(te.fcc(2), list(FCC2_LARGE), label_exchange=True) == FCC2_LARGE


@pytest.mark.parametrize("name,k,nmax", [("fcc", 2, 10), ("bcc", 3, 7), ("sc", 2, 9), ("sc", 4, 5),
                                         ("hcp", 2, 6), ("hcp", 3, 4)])
@pytest.mark.parametrize("label_exchange,include_incomplete", list(itertools.product((False, True), repeat=2)))
def test_matches_generator(name, k, nmax, label_exchange, include_incomplete):
    P = getattr(te, name)(k)
    _agree(P, range(1, nmax + 1), label_exchange=label_exchange, include_incomplete=include_incomplete)


@pytest.mark.parametrize("P,n", [(te.fcc(2), 8), (te.hcp(2), 5), (te.fcc(3), 6)])
@pytest.mark.parametrize("label_exchange", [False, True])
def test_fixed_compositions(P, n, label_exchange):
    N = n * P.d
    for cut in itertools.combinations(range(N + 1), P.k - 1):
        c = tuple(b - a for a, b in zip((0,) + cut, cut + (N,)))
        _agree(P, n, composition=c, label_exchange=label_exchange)


@pytest.mark.parametrize("allowed", [[{"Li", "Vac"}, {"Co", "Vac"}], [{"Li"}, {"Vac", "Co"}]])
@pytest.mark.parametrize("label_exchange", [False, True])
def test_site_restrictions(allowed, label_exchange):
    P = te.hcp(["Li", "Vac", "Co"], allowed=allowed)
    _agree(P, range(1, 6), label_exchange=label_exchange)
    _agree(P, range(1, 6), label_exchange=label_exchange, include_incomplete=True)
