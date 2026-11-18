"""Independent brute-force reference.

Shares no code with the word machinery: sites are lattice points folded into
the supercell by coordinates, symmetry operations are applied to Cartesian
positions, superlattices are checked with the pairwise test of Hart & Forcade,
and every allowed labeling is canonicalized over the full group.
"""
import itertools
from collections import defaultdict

import numpy as np

from necklace_enum.lattice import Quotient, hnfs


def is_int(X):
    return np.allclose(X, np.round(X), atol=1e-8)


def label_group(P, label_exchange):
    if not label_exchange:
        return [tuple(range(P.k))]
    return [p for p in itertools.permutations(range(P.k))
            if all(frozenset(p[l] for l in s) == s for s in P.allowed)]


def check_superlattices(P, n, reps):
    def equiv(H1, H2):
        H2i = np.linalg.inv(H2)
        return any(is_int(H2i @ M @ H1) for M in P.rots)
    for i, j in itertools.combinations(range(len(reps)), 2):
        assert not equiv(reps[i], reps[j]), "duplicate superlattice representatives"
    for H in hnfs(n):
        assert any(equiv(H, R) for R in reps), "an HNF is not covered"


class Cell:
    def __init__(self, P, H):
        self.P, self.H = P, H
        self.Hi = np.linalg.inv(H)
        keys, zs = {}, []
        r = int(np.abs(H).max()) + 1
        for z in itertools.product(range(-r, r + 1), repeat=3):
            k_ = self.key(np.array(z))
            if k_ not in keys:
                keys[k_] = len(zs)
                zs.append(np.array(z))
        assert len(zs) == round(abs(np.linalg.det(H)))
        self.keys, self.zs, self.n = keys, zs, len(zs)

    def key(self, z):
        f = self.Hi @ z
        f = f - np.floor(f + 1e-9)
        return tuple(np.round(f, 6) % 1.0)

    def site(self, i, z):
        return i * self.n + self.keys[self.key(z)]

    def group(self):
        P = self.P
        trans = [[self.site(i, z + t) for i in range(P.d) for z in self.zs] for t in self.zs]
        ops = []
        for M, tau, pi, v in P.ops:
            if not is_int(self.Hi @ M @ self.H):
                continue
            R = P.A @ M @ P.Ainv
            p = []
            for i in range(P.d):
                for z in self.zs:
                    x = R @ (P.A @ z + P.basis[i]) + tau
                    for j in range(P.d):
                        f = P.Ainv @ (x - P.basis[j])
                        if is_int(f):
                            p.append(self.site(j, np.round(f).astype(int)))
                            break
                    else:
                        raise AssertionError("operation does not map sites to sites")
            ops.append(p)
        N = self.n * P.d
        zero = self.keys[self.key(np.zeros(3))]
        nontrivial = [t for idx, t in enumerate(trans) if idx != zero]
        return [[o[t[s]] for s in range(N)] for o in ops for t in trans], nontrivial


def canon(lab, group, LG):
    best = None
    for g in group:
        for e in LG:
            img = [0] * len(lab)
            for s in range(len(lab)):
                img[g[s]] = e[lab[s]]
            img = tuple(img)
            if best is None or img < best:
                best = img
    return best


def brute(P, n, reps, label_exchange, conc):
    LG = label_group(P, label_exchange)
    N = n * P.d
    out = {}
    for H in reps:
        cell = Cell(P, H)
        group, nontrivial = cell.group()
        choices = [sorted(P.allowed[s // n]) for s in range(N)]
        found = set()
        for lab in itertools.product(*choices):
            cnt = [lab.count(l) for l in range(P.k)]
            if conc is not None and cnt != list(conc):
                continue
            if conc is None and min(cnt) == 0:
                continue
            if any(tuple(lab[t[s]] for s in range(N)) == lab for t in nontrivial):
                continue
            found.add(canon(lab, group, LG))
        out[tuple(H.ravel())] = (cell, group, found)
    return out


def to_cell_labeling(P, cell, H, lab):
    q = Quotient(H)
    new = [0] * (q.n * P.d)
    for i in range(P.d):
        for g in range(q.n):
            new[cell.site(i, q.Linv @ q.elems[g])] = lab[i * q.n + g]
    return new
