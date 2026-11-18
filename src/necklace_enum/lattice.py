"""Superlattices, Smith normal forms and quotient groups.

Conventions follow Hart & Forcade, Phys. Rev. B 77, 224115 (2008). Everything
is in lattice coordinates: a parent lattice is a 3x3 matrix A whose *columns*
are the primitive vectors, a superlattice is A @ H with H an integer matrix in
lower-triangular Hermite normal form, and point operations are integer
matrices R with R^T G R = G for the metric G = A^T A.
"""
from __future__ import annotations

import itertools

import numpy as np
from numba import njit

def point_group(A: np.ndarray, tol: float = 1e-8) -> list[np.ndarray]:
    """All integer matrices R (lattice coords) that preserve the metric."""
    G = A.T @ A
    lengths = np.diag(G)
    # candidate images of each basis vector: lattice vectors of equal length
    rng = range(-3, 4)
    pts = np.array(list(itertools.product(rng, rng, rng)))
    norms = np.einsum("ij,jk,ik->i", pts, G, pts)
    cands = [pts[np.abs(norms - L) < tol] for L in lengths]
    ops = []
    for v0 in cands[0]:
        for v1 in cands[1]:
            if abs(v0 @ G @ v1 - G[0, 1]) > tol:
                continue
            for v2 in cands[2]:
                R = np.array([v0, v1, v2]).T
                if np.allclose(R.T @ G @ R, G, atol=tol):
                    ops.append(R.astype(np.int64))
    return ops


# ---------------------------------------------------------------- HNF / SNF
def hnfs(n: int) -> list[np.ndarray]:
    """All lower-triangular HNFs of determinant n (paper eq. 1)."""
    out = []
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
                        out.append(np.array([[a, 0, 0], [b, c, 0], [d, e, f]], dtype=np.int64))
    return out


def to_hnf(M: np.ndarray) -> np.ndarray:
    """Lower-triangular HNF of M by integer column operations."""
    M = M.astype(np.int64).copy()
    for i in range(3):
        # zero M[i, j] for j > i using column ops on columns i..2
        for j in range(i + 1, 3):
            while M[i, j] != 0:
                q = M[i, i] // M[i, j]
                M[:, i] -= q * M[:, j]
                M[:, [i, j]] = M[:, [j, i]]
        if M[i, i] < 0:
            M[:, i] *= -1
    for i in range(1, 3):
        for j in range(i):
            q = M[i, j] // M[i, i]  # floor: gives 0 <= M[i,j] < M[i,i]
            M[:, j] -= q * M[:, i]
    return M


def snf(H: np.ndarray):
    """Return (D, L, R) with D = L @ H @ R diagonal, d1 | d2 | d3, all > 0,
    L and R unimodular integer matrices."""
    S = H.astype(np.int64).copy()
    L = np.eye(3, dtype=np.int64)
    R = np.eye(3, dtype=np.int64)
    for t in range(3):
        while True:
            # move smallest nonzero |entry| of the trailing block to (t,t)
            sub = np.abs(S[t:, t:])
            if not sub.any():
                break
            sub = np.where(sub == 0, np.iinfo(np.int64).max, sub)
            i, j = np.unravel_index(np.argmin(sub), sub.shape)
            i += t; j += t
            S[[t, i]] = S[[i, t]]; L[[t, i]] = L[[i, t]]
            S[:, [t, j]] = S[:, [j, t]]; R[:, [t, j]] = R[:, [j, t]]
            done = True
            for i in range(t + 1, 3):  # clear column t with row ops
                q = S[i, t] // S[t, t]
                S[i] -= q * S[t]; L[i] -= q * L[t]
                if S[i, t]:
                    done = False
            for j in range(t + 1, 3):  # clear row t with column ops
                q = S[t, j] // S[t, t]
                S[:, j] -= q * S[:, t]; R[:, j] -= q * R[:, t]
                if S[t, j]:
                    done = False
            if not done:
                continue
            # divisibility: pivot must divide the rest of the block
            bad = [(i, j) for i in range(t + 1, 3) for j in range(t + 1, 3)
                   if S[i, j] % S[t, t]]
            if bad:
                i, _ = bad[0]
                S[t] += S[i]; L[t] += L[i]
                continue
            break
        if S[t, t] < 0:
            S[t] *= -1; L[t] *= -1
    assert np.array_equal(L @ H @ R, S)
    assert abs(round(np.linalg.det(L))) == 1 and abs(round(np.linalg.det(R))) == 1
    d = np.diag(S)
    assert d[0] and d[1] % d[0] == 0 and d[2] % d[1] == 0
    return S, L, R


def stabilizer(H, rots):
    """Operations R (lattice coordinates) that map the superlattice H onto itself."""
    Hi = np.linalg.inv(H)
    out = []
    for R in rots:
        X = Hi @ R @ H
        if np.allclose(X, np.round(X), atol=1e-8):
            out.append(R)
    return out


# ---------------------------------------------------------------- quotient group
class Quotient:
    """Quotient group L/S = Z_d1 + Z_d2 + Z_d3 for superlattice H, with the
    permutations of its n elements induced by translations and rotations."""

    def __init__(self, H: np.ndarray):
        self.H = H
        self.D, self.L, self.R = snf(H)
        self.d = tuple(int(x) for x in np.diag(self.D))
        self.n = int(np.prod(self.d))
        self.elems = np.array(list(itertools.product(*(range(x) for x in self.d))),
                              dtype=np.int64)  # mixed-radix order
        self.Linv = np.round(np.linalg.inv(self.L)).astype(np.int64)

    def index(self, g: np.ndarray) -> np.ndarray:
        g = np.mod(g, self.d)
        return (g[..., 0] * self.d[1] + g[..., 1]) * self.d[2] + g[..., 2]

    def translation_perms(self) -> np.ndarray:
        """perm[t][s] = index of (elem_s + elem_t)."""
        E = self.elems
        return self.index(E[None, :, :] + E[:, None, :]).astype(np.int64)

    def rotation_perm(self, R: np.ndarray) -> np.ndarray:
        """g -> L R L^{-1} g  (mod D)."""
        M = self.L @ R @ self.Linv
        p = self.index(self.elems @ M.T).astype(np.int64)
        assert len(set(p.tolist())) == self.n, "rotation does not fix superlattice"
        return p


    def column_major(self):
        """phi[g] = column-major position (column g3, then track g1*d2 + g2),
        the number of tracks a = d1*d2 and the column count m = d3."""
        d1, d2, m = self.d
        a = d1 * d2
        E = self.elems
        return (E[:, 2] * a + E[:, 0] * d2 + E[:, 1]).astype(np.int64), a, m


@njit(cache=True)
def _hnf_key(M, base):
    """Lower-triangular HNF of 3x3 integer M (column ops), packed into one int."""
    M = M.copy()
    for i in range(3):
        for j in range(i + 1, 3):
            while M[i, j] != 0:
                q = M[i, i] // M[i, j]
                for r in range(3):
                    M[r, i] -= q * M[r, j]
                for r in range(3):
                    M[r, i], M[r, j] = M[r, j], M[r, i]
        if M[i, i] < 0:
            for r in range(3):
                M[r, i] = -M[r, i]
    for i in range(1, 3):
        for j in range(i):
            q = M[i, j] // M[i, i]
            for r in range(3):
                M[r, j] -= q * M[r, i]
    key = 0
    for (i, j) in ((0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)):
        key = key * base + M[i, j]
    return key


@njit(cache=True)
def _canonical_keys(Hs, Rs, base):
    out = np.empty(Hs.shape[0], dtype=np.int64)
    P = np.empty((3, 3), dtype=np.int64)
    for h in range(Hs.shape[0]):
        best = -1
        for r in range(Rs.shape[0]):
            for i in range(3):
                for j in range(3):
                    s = 0
                    for l in range(3):
                        s += Rs[r, i, l] * Hs[h, l, j]
                    P[i, j] = s
            k = _hnf_key(P, base)
            if best < 0 or k < best:
                best = k
        out[h] = best
    return out


def unique_superlattices(n, rots):
    """Symmetry-distinct HNFs of index n. Each class is represented by the
    lexicographically smallest HNF of R @ H over the operations R."""
    Hs = np.array(hnfs(n), dtype=np.int64)
    Rs = np.array(rots, dtype=np.int64)
    base = n + 1
    keys = _canonical_keys(Hs, Rs, base)
    _, first = np.unique(keys, return_index=True)
    reps = []
    for key in keys[np.sort(first)]:
        v = []
        for _ in range(6):
            v.append(int(key % base)); key //= base
        a, b, c, d, e, f = reversed(v)
        reps.append(np.array([[a, 0, 0], [b, c, 0], [d, e, f]], dtype=np.int64))
    return reps
