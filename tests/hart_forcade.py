"""Reference implementation of the Hart & Forcade (2008) algorithm, for comparison.

Follows Phys. Rev. B 77, 224115 (2008), Sec. III, for Bravais parents:
  3a  enumerate all k^n labelings of each Smith normal form, indexed by a
      perfect hash (the base-k rank of the labeling)
  3b  drop labelings in which some species is absent
  3c,d  mark every image under translations and label exchange as a duplicate
  3e  drop labelings that are periodic in a smaller cell
  4   for each superlattice, drop duplicates under the point operations that fix it
      (images are mapped back to their step-3 representative before marking)

The kernels are compiled with Numba, like necklace-enum, so the timings compare
algorithms rather than languages. Memory is 8 k^n bytes per Smith normal form.
"""
from collections import defaultdict

import numpy as np
from numba import njit

from necklace_enum.lattice import Quotient, stabilizer, unique_superlattices

INCOMPLETE, SUPERPER = -2, -3


@njit(cache=True)
def _rank(lab, k):
    r = 0
    for s in range(lab.shape[0]):
        r = r * k + lab[s]
    return r


@njit(cache=True)
def _unrank(x, k, n, out):
    for s in range(n - 1, -1, -1):
        out[s] = x % k
        x //= k


@njit(cache=True)
def _snf_stage(k, n, tperms, lperms):
    """Steps 3a-3e. rep[x] == x marks a survivor; otherwise rep[x] is its representative."""
    N = k ** n
    rep = np.full(N, -1, dtype=np.int64)
    lab = np.empty(n, dtype=np.int64)
    img = np.empty(n, dtype=np.int64)
    seen = np.empty(k, dtype=np.bool_)
    for x in range(N):
        if rep[x] != -1:
            continue
        _unrank(x, k, n, lab)
        seen[:] = False
        for s in range(n):
            seen[lab[s]] = True
        if not seen.all():
            rep[x] = INCOMPLETE
            continue
        rep[x] = x
        superper = False
        for t in range(tperms.shape[0]):
            for e in range(lperms.shape[0]):
                for s in range(n):
                    img[tperms[t, s]] = lperms[e, lab[s]]
                y = _rank(img, k)
                if y == x and t != 0 and e == 0:
                    superper = True
                if rep[y] == -1:
                    rep[y] = x
        if superper:
            rep[x] = SUPERPER
    return rep


@njit(cache=True)
def _hnf_stage(k, n, survivors, rep, rperms):
    """Step 4 for one superlattice: survivors that are unique under its point operations."""
    dup = np.zeros(rep.shape[0], dtype=np.bool_)
    lab = np.empty(n, dtype=np.int64)
    img = np.empty(n, dtype=np.int64)
    m = 0
    for i in range(survivors.shape[0]):
        x = survivors[i]
        if dup[x]:
            continue
        m += 1
        _unrank(x, k, n, lab)
        for r in range(rperms.shape[0]):
            for s in range(n):
                img[rperms[r, s]] = lab[s]
            y = _rank(img, k)
            y = rep[y] if rep[y] >= 0 else y
            if y != x:
                dup[y] = True
    return m


def count_structures(parent, n, label_exchange=True):
    """Number of derivative structures of index n (Bravais parent, all compositions,
    every species present), by the Hart-Forcade algorithm."""
    if parent.d != 1:
        raise ValueError("reference implementation covers Bravais parents only")
    k = parent.k
    import itertools
    lp = np.array(list(itertools.permutations(range(k))) if label_exchange else [tuple(range(k))],
                  dtype=np.int64)
    by_snf = defaultdict(list)
    for H in unique_superlattices(n, parent.rots):
        q = Quotient(H)
        by_snf[q.d].append(q)
    total = 0
    for qs in by_snf.values():
        rep = _snf_stage(k, n, qs[0].translation_perms(), lp)
        survivors = np.nonzero(rep == np.arange(rep.size))[0]
        for q in qs:
            rperms = np.array([q.rotation_perm(R) for R in stabilizer(q.H, parent.rots)], dtype=np.int64)
            total += _hnf_stage(k, n, survivors, rep, rperms)
    return total
