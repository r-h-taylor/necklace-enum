"""Structure-by-structure comparison with the Hart-Forcade reference implementation.

For every case, both codes list the structures of one cell size, superlattice by
superlattice. Each labeling is reduced to a canonical form, the smallest labeling in
its orbit under the superlattice's translations, its point operations and (with label
exchange) the permutations of the species. The two lists must be equal for every
superlattice, with no duplicates.
"""
import itertools

import numpy as np
import pytest

import necklace_enum as te
from necklace_enum.lattice import Quotient, stabilizer

import hart_forcade as hf

SQ3 = np.sqrt(3) / 2
PARENTS = {
    "fcc": te.fcc,
    "bcc": te.bcc,
    "sc": te.sc,
    "hex": lambda k: te.Parent([[1, 0, 0], [0.5, SQ3, 0], [0, 0, 1.6]], species=k),
    "tet": lambda k: te.Parent([[1, 0, 0], [0, 1, 0], [0, 0, 1.3]], species=k),
}

# (parent, k, n, composition or None); each runs with and without label exchange
SIZES = {
    ("fcc", 2): range(1, 11), ("fcc", 3): range(1, 7), ("fcc", 4): range(1, 6),
    ("bcc", 2): range(1, 9), ("bcc", 3): range(1, 5),
    ("sc", 2): range(1, 9), ("sc", 3): range(1, 5),
    ("hex", 2): range(1, 7), ("hex", 3): range(1, 5),
    ("tet", 2): range(1, 7), ("tet", 3): range(1, 5),
}
FIXED = [("fcc", 2, 8, (6, 2)), ("fcc", 2, 8, (4, 4)), ("fcc", 2, 10, (7, 3)), ("fcc", 3, 6, (3, 2, 1)),
         ("fcc", 4, 6, (2, 2, 1, 1)), ("bcc", 2, 8, (5, 3)), ("sc", 2, 8, (6, 2)), ("sc", 3, 6, (2, 2, 2)),
         ("hex", 2, 6, (4, 2)), ("hex", 3, 4, (2, 1, 1)), ("tet", 2, 6, (3, 3)), ("tet", 3, 4, (2, 1, 1))]
CASES = ([(p, k, n, None, le) for (p, k), ns in SIZES.items() for n in ns if n >= k for le in (False, True)]
         + [(p, k, n, c, le) for p, k, n, c in FIXED for le in (False, True)])
SLOW = {("fcc", 2, 10), ("fcc", 2, 9)}


def canonical(H, lab, rots, lperms):
    """Smallest labeling in the orbit of lab on superlattice H."""
    q = Quotient(H)
    lab = np.asarray(lab)
    T = q.translation_perms()
    best = None
    for R in stabilizer(H, rots):
        rp = q.rotation_perm(R)
        r_img = np.empty_like(lab)
        r_img[rp] = lab
        imgs = np.empty((T.shape[0], lab.size), dtype=lab.dtype)
        for t in range(T.shape[0]):
            imgs[t, T[t]] = r_img
        for e in lperms:
            for row in e[imgs]:
                key = tuple(int(v) for v in row)
                if best is None or key < best:
                    best = key
    return best


def lists(structures, rots, lperms):
    out = {}
    for H, lab in structures:
        out.setdefault(tuple(np.asarray(H).ravel()), []).append(canonical(H, lab, rots, lperms))
    return out


@pytest.mark.parametrize("parent,k,n,comp,le", [
    pytest.param(*c, marks=pytest.mark.slow) if c[:3] in SLOW else c for c in CASES])
def test_lists_agree(parent, k, n, comp, le):
    P = PARENTS[parent](k)
    lperms = np.array(list(itertools.permutations(range(k))) if le else [tuple(range(k))])
    ref = lists(hf.enumerate_structures(P, n, label_exchange=le, composition=comp), P.rots, lperms)
    got = lists(((s.hnf, s.labels) for s in te.enumerate_structures(P, n, composition=comp, label_exchange=le)),
                P.rots, lperms)
    for h in set(ref) | set(got):
        mine, theirs = got.get(h, []), ref.get(h, [])
        assert len(mine) == len(set(mine)), f"duplicate structures on {h}"
        assert len(theirs) == len(set(theirs)), f"reference duplicates on {h}"
        assert set(mine) == set(theirs), f"lists differ on {h}"


@pytest.mark.parametrize("k,nmax", [(2, 12), (3, 8), (4, 7)])
def test_counts_agree(k, nmax):
    P = te.fcc(k)
    for n in range(2, nmax + 1):
        assert hf.count_structures(P, n) == te.count_structures(P, n, label_exchange=True)[n]

