"""Degeneracies computed inside the generator.

Each structure's degeneracy (enumlib's labdegn) comes from counting the symmetry comparisons that
end equal while the word is tested for canonicity: orbit size = group order / stabilizer order.
It is checked against (1) an independent post-pass that applies every site permutation of the
superlattice's group to the labeling, and (2) a direct count of every labeling of every superlattice:
summed over all structures of one size, the total degeneracies must count each once.
"""
import itertools

import numpy as np
import pytest

import necklace_enum as ne
from necklace_enum.core import _label_group
from necklace_enum.enumlib_format import _HNFInfo
from necklace_enum.lattice import Quotient
from necklace_enum.parent import Parent

ROCKSALT = Parent([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]], [[0, 0, 0], [.5, .5, .5]], species=3,
                  allowed=[[0, 1], [2]])
DIAMOND = Parent([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]], [[0, 0, 0], [.25, .25, .25]], species=2)

CASES = [
    ("fcc", ne.fcc(2), range(1, 10), {}),
    ("fcc-lex", ne.fcc(2), range(1, 10), {"label_exchange": True}),
    ("fcc-incomplete", ne.fcc(2), range(1, 8), {"include_incomplete": True}),
    ("fcc3", ne.fcc(3), range(1, 7), {}),
    ("fcc3-lex", ne.fcc(3), range(1, 7), {"label_exchange": True}),
    ("bcc", ne.bcc(2), range(1, 9), {}),
    ("sc3-lex", ne.sc(3), range(1, 6), {"label_exchange": True}),
    ("hcp", ne.hcp(2), range(1, 6), {}),
    ("hcp-lex", ne.hcp(2), range(1, 6), {"label_exchange": True}),
    ("fcc-6:2", ne.fcc(2), [8], {"composition": (6, 2)}),
    ("fcc-4:4-lex", ne.fcc(2), [8], {"composition": (4, 4), "label_exchange": True}),
    ("fcc-3:3:2", ne.fcc(3), [8], {"composition": (3, 3, 2)}),
    ("hcp-5:3", ne.hcp(2), [4], {"composition": (5, 3)}),
    ("rocksalt-restricted", ROCKSALT, range(1, 6), {}),
    ("diamond-lex", DIAMOND, range(1, 5), {"label_exchange": True}),
]


@pytest.mark.parametrize("name,P,sizes,kw", CASES, ids=[c[0] for c in CASES])
def test_generator_degeneracy_matches_post_pass(name, P, sizes, kw):
    LG = _label_group(P, kw.get("label_exchange", False))
    for n in sizes:
        info = {}
        for H, labs, degs in ne.enumerate_labelings(P, n, degeneracy=True, **kw):
            key = tuple(np.asarray(H).ravel())
            h = info.get(key) or info.setdefault(key, _HNFInfo(P, np.asarray(H), 0))
            assert (h.labdegn(labs.astype(np.int64), LG) == degs).all()


@pytest.mark.parametrize("name,P,sizes,kw", CASES, ids=[c[0] for c in CASES])
def test_structure_objects_and_counting_agree(name, P, sizes, kw):
    for n in sizes:
        st = {}
        ne.count_structures(P, n, stats=st, **kw)
        structs = list(ne.enumerate_structures(P, n, **kw))
        assert all(s.total_degeneracy == s.degeneracy * s.superlattice_degeneracy for s in structs)
        assert all(s.superlattice_degeneracy == ne.superlattice_degeneracy(P, s.hnf) for s in structs)
        assert st.get("weight", 0) == sum(s.total_degeneracy for s in structs)


def _all_hnfs(n):
    for a in range(1, n + 1):
        if n % a:
            continue
        for c in range(1, n // a + 1):
            if (n // a) % c:
                continue
            f = n // (a * c)
            for b in range(c):
                for d in range(f):
                    for e in range(f):
                        yield np.array([[a, 0, 0], [b, c, 0], [d, e, f]])


def _every_labeling(P, n, comp):
    """Labelings of every superlattice of index n, not periodic in a smaller cell."""
    k, d, tot = P.k, P.d, 0
    for H in _all_hnfs(n):
        tp = np.asarray(Quotient(H).translation_perms())
        full = [np.concatenate([i * n + tp[t] for i in range(d)]) for t in range(1, n)]
        for lab in itertools.product(range(k), repeat=n * d):
            lab = np.array(lab)
            cnt = np.bincount(lab, minlength=k)
            if (comp is None and cnt.min() == 0) or (comp is not None and tuple(cnt) != comp):
                continue
            tot += not any((lab[f] == lab).all() for f in full)
    return tot


@pytest.mark.parametrize("P,n,comp,lex", [(ne.fcc(2), 8, None, False), (ne.fcc(2), 8, None, True),
                                          (ne.fcc(3), 5, None, False), (ne.hcp(2), 4, None, False),
                                          (ne.hcp(2), 4, (5, 3), False), (ne.fcc(2), 8, (5, 3), False)])
def test_total_weight_counts_every_labeling(P, n, comp, lex):
    st = {}
    ne.count_structures(P, n, composition=comp, label_exchange=lex, stats=st)
    assert st["weight"] == _every_labeling(P, n, comp)
