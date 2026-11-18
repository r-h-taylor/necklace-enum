"""Parent lattices: Bravais lattices and multilattices with site restrictions."""
from __future__ import annotations

import numpy as np

from .lattice import point_group


class Parent:
    """A parent crystal on which derivative structures are built.

    Parameters
    ----------
    lattice : (3, 3) array
        Primitive lattice vectors as *rows*, in Cartesian coordinates.
    basis : (d, 3) array
        Cartesian positions of the d sites in the primitive cell.
    species : sequence of str or int
        The atom types (e.g. ``["Cu", "Au"]``). An integer k means k unnamed
        types labelled 0..k-1.
    allowed : sequence of sequences, optional
        For each basis site, the species allowed on it (names or indices).
        Default: every species on every site.
    """

    def __init__(self, lattice, basis=((0.0, 0.0, 0.0),), species=2, allowed=None, eps=1e-6):
        if isinstance(species, int):
            self.species = list(range(species))
        else:
            self.species = list(species)
        self.k = len(self.species)
        if self.k < 2:
            raise ValueError("need at least two species")
        self.A = np.array(lattice, dtype=float).T          # columns = vectors
        self.Ainv = np.linalg.inv(self.A)
        self.basis = [np.array(x, dtype=float) for x in basis]
        self.d = len(self.basis)
        if allowed is None:
            allowed = [range(self.k)] * self.d
        if len(allowed) != self.d:
            raise ValueError("allowed must have one entry per basis site")
        self.allowed = [frozenset(self._index(s) for s in a) for a in allowed]
        self.eps = eps
        self.ops = self._space_group()
        self.rots, seen = [], set()
        for M, *_ in self.ops:
            if M.tobytes() not in seen:
                seen.add(M.tobytes())
                self.rots.append(M)

    def _index(self, s):
        if isinstance(s, (int, np.integer)) and s not in self.species:
            if not 0 <= s < self.k:
                raise ValueError(f"species index {s} out of range")
            return int(s)
        return self.species.index(s)

    def _space_group(self):
        """Operations (M, tau, pi, v): rotation M in lattice coordinates, Cartesian
        translation tau, site-kind permutation pi and lattice offsets v[i] with
        R d_i + tau = d_pi(i) + A v_i. Operations must preserve site restrictions."""
        ops = []
        for M in point_group(self.A):
            R = self.A @ M @ self.Ainv
            RD = [R @ x for x in self.basis]
            for j in range(self.d):
                tau = self.basis[j] - RD[0]
                pi, v, ok = [], [], True
                for i in range(self.d):
                    hit = None
                    for jj in range(self.d):
                        f = self.Ainv @ (RD[i] + tau - self.basis[jj])
                        if np.allclose(f, np.round(f), atol=self.eps):
                            hit = (jj, np.round(f).astype(np.int64))
                            break
                    if hit is None or self.allowed[hit[0]] != self.allowed[i]:
                        ok = False
                        break
                    pi.append(hit[0])
                    v.append(hit[1])
                if ok:
                    ops.append((M, tau, tuple(pi), np.array(v)))
        uniq, keys = [], set()
        for M, tau, pi, v in ops:
            key = (M.tobytes(), pi, v.tobytes())
            if key not in keys:
                keys.add(key)
                uniq.append((M, tau, pi, v))
        return uniq

    def __repr__(self):
        return f"Parent(d={self.d}, species={self.species}, {len(self.ops)} space-group operations)"


# ---------------------------------------------------------------- common parents
def fcc(species=2, a=1.0):
    return Parent(a * np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]), species=species)


def bcc(species=2, a=1.0):
    return Parent(a * np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]), species=species)


def sc(species=2, a=1.0):
    return Parent(a * np.eye(3), species=species)


def hcp(species=2, a=1.0, c=None, allowed=None):
    c = np.sqrt(8.0 / 3.0) * a if c is None else c
    lat = [[a, 0, 0], [0.5 * a, np.sqrt(3) / 2 * a, 0], [0, 0, c]]
    basis = [[0, 0, 0], [0.5 * a, np.sqrt(3) / 6 * a, c / 2]]
    return Parent(lat, basis, species=species, allowed=allowed)
