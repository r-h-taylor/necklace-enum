"""Reference implementation of the Hart & Forcade (2008) algorithm, for comparison.

Follows Phys. Rev. B 77, 224115 (2008), Sec. III, for Bravais parents:
  3a  enumerate all k^n labelings of each Smith normal form, indexed by a
      perfect hash (the base-k rank of the labeling)
  3b  drop labelings in which some species is absent, and, at a fixed
      composition, labelings with any other composition
  3c,d  mark every image under translations and label exchange as a duplicate
  3e  drop labelings that are periodic in a smaller cell
  4   for each superlattice, drop duplicates under the point operations that fix it
      (images are mapped back to their step-3 representative before marking;
      a literal reading of step 4, which compares r.x only with the step-3
      survivors, overcounts from n = 7 on for binary fcc)

It shares the superlattice list and the Smith normal forms (necklace_enum.lattice)
with the generator, but none of the word machinery. The kernels are compiled with
Numba. Memory is 8 k^n bytes per Smith normal form, so it is for small cells only.
"""
import itertools
from collections import defaultdict

import numpy as np
from numba import njit

from necklace_enum.lattice import Quotient, stabilizer, unique_superlattices

INCOMPLETE, SUPERPER, OTHER = -2, -3, -4


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
def _snf_stage(k, n, tperms, lperms, conc):
    """Steps 3a-3e. rep[x] == x marks a survivor; otherwise rep[x] is its representative
    or a negative code. conc[0] < 0 means all compositions."""
    N = k ** n
    rep = np.full(N, -1, dtype=np.int64)
    lab = np.empty(n, dtype=np.int64)
    img = np.empty(n, dtype=np.int64)
    cnt = np.zeros(k, dtype=np.int64)
    for x in range(N):
        if rep[x] != -1:
            continue
        _unrank(x, k, n, lab)
        cnt[:] = 0
        for s in range(n):
            cnt[lab[s]] += 1
        if conc[0] >= 0:
            same = True
            for l in range(k):
                if cnt[l] != conc[l]:
                    same = False
            if not same:
                rep[x] = OTHER
                continue
        incomplete = False
        for l in range(k):
            if cnt[l] == 0:
                incomplete = True
        if incomplete:
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
                if rep[y] == -1 or rep[y] == OTHER:
                    rep[y] = x
        if superper:
            rep[x] = SUPERPER
    return rep


@njit(cache=True)
def _hnf_stage(k, n, survivors, rep, rperms):
    """Step 4 for one superlattice: the survivors that are unique under its point operations."""
    dup = np.zeros(rep.shape[0], dtype=np.bool_)
    keep = np.zeros(survivors.shape[0], dtype=np.bool_)
    lab = np.empty(n, dtype=np.int64)
    img = np.empty(n, dtype=np.int64)
    for i in range(survivors.shape[0]):
        x = survivors[i]
        if dup[x]:
            continue
        keep[i] = True
        _unrank(x, k, n, lab)
        for r in range(rperms.shape[0]):
            for s in range(n):
                img[rperms[r, s]] = lab[s]
            y = _rank(img, k)
            y = rep[y] if rep[y] >= 0 else y
            if y != x:
                dup[y] = True
    return survivors[keep]


def _run(parent, n, label_exchange, composition):
    if parent.d != 1:
        raise ValueError("reference implementation covers Bravais parents only")
    k = parent.k
    lp = np.array(list(itertools.permutations(range(k))) if label_exchange else [tuple(range(k))],
                  dtype=np.int64)
    conc = np.full(k, -1, dtype=np.int64) if composition is None else np.array(composition, dtype=np.int64)
    if composition is not None and conc.sum() != n:
        raise ValueError("composition must add up to n")
    by_snf = defaultdict(list)
    for H in unique_superlattices(n, parent.rots):
        q = Quotient(H)
        by_snf[q.d].append(q)
    lab = np.empty(n, dtype=np.int64)
    for qs in by_snf.values():
        rep = _snf_stage(k, n, qs[0].translation_perms(), lp, conc)
        survivors = np.nonzero(rep == np.arange(rep.size))[0]
        for q in qs:
            rperms = np.array([q.rotation_perm(R) for R in stabilizer(q.H, parent.rots)], dtype=np.int64)
            kept = _hnf_stage(k, n, survivors, rep, rperms)
            out = []
            for x in kept:
                _unrank(x, k, n, lab)
                out.append(tuple(int(v) for v in lab))
            yield q.H, out


def count_structures(parent, n, label_exchange=True, composition=None):
    """Number of derivative structures of index n (Bravais parent, every species present,
    all compositions or one composition given as atom counts), by the Hart-Forcade algorithm."""
    return sum(len(labs) for _, labs in _run(parent, n, label_exchange, composition))


def enumerate_structures(parent, n, label_exchange=True, composition=None):
    """(hnf, labeling) for every structure of index n, labeling indexed like
    necklace_enum.Structure.labels (quotient-group elements in mixed-radix order)."""
    for H, labs in _run(parent, n, label_exchange, composition):
        for lab in labs:
            yield H, lab
