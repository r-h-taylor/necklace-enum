"""Enumeration of derivative structures as canonical words."""
from __future__ import annotations

import itertools
from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from .kernel import product_tables, reversal_maps, split_reversal, torus_kernel, torus_kernel_resume
from .lattice import Quotient, unique_superlattices
from .parent import Parent

LIST_BATCH = 1 << 16          # words per batch when listing; memory stays fixed at about this size


# ---------------------------------------------------------------- results
@dataclass(frozen=True)
class Structure:
    """One derivative structure.

    ``hnf`` is the superlattice in Hermite normal form (columns are supercell
    vectors in lattice coordinates), and ``labels[i * n + g]`` is the species
    index on basis site i at quotient-group element g (mixed-radix order of the
    Smith normal form).
    """
    parent: Parent
    hnf: np.ndarray
    labels: tuple

    @property
    def size(self):
        return int(round(abs(np.linalg.det(self.hnf))))

    @property
    def cell(self):
        """Supercell vectors as rows, Cartesian."""
        return (self.parent.A @ self.hnf).T

    def sites(self):
        """Cartesian positions and species of every atom, wrapped into the supercell."""
        q = Quotient(self.hnf)
        A, n = self.parent.A, q.n
        z = q.elems @ q.Linv.T                              # lattice point with g(z) = g
        pos, spc = [], []
        for i, b in enumerate(self.parent.basis):
            for g in range(n):
                pos.append(A @ z[g] + b)
                spc.append(self.parent.species[self.labels[i * n + g]])
        pos = np.array(pos)
        frac = np.linalg.solve(self.cell.T, pos.T).T
        frac -= np.floor(frac + 1e-9)
        return frac @ self.cell, spc

    def to_ase(self):
        from ase import Atoms
        pos, spc = self.sites()
        return Atoms(symbols=[str(s) for s in spc], positions=pos, cell=self.cell, pbc=True)

    def to_pymatgen(self):
        from pymatgen.core import Lattice, Structure as PmgStructure
        pos, spc = self.sites()
        return PmgStructure(Lattice(self.cell), [str(s) for s in spc], pos, coords_are_cartesian=True)


# ---------------------------------------------------------------- setup helpers
def _site_perm(P: Parent, q: Quotient, op):
    M, tau, pi, v = op
    n, d = q.n, P.d
    T = q.L @ M @ q.Linv
    perm = np.empty(n * d, dtype=np.int64)
    for i in range(d):
        img = q.index(q.elems @ T.T + q.L @ v[i])
        perm[i * n:(i + 1) * n] = pi[i] * n + img
    return perm


def _stabilizer_ops(P: Parent, H):
    Hi = np.linalg.inv(H)
    out = []
    for op in P.ops:
        X = Hi @ op[0] @ H
        if np.allclose(X, np.round(X), atol=1e-8):
            out.append(op)
    return out


def _column_major(q: Quotient, d):
    """phi[site] = word position: column g3, then site kind, then track."""
    phi1, a, m = q.column_major()
    n = q.n
    phi = np.empty(n * d, dtype=np.int64)
    g3, tr = phi1 // a, phi1 % a
    for i in range(d):
        phi[i * n:(i + 1) * n] = g3 * (d * a) + i * a + tr
    return phi, a * d, m


def _label_group(P: Parent, label_exchange):
    if not label_exchange:
        return [np.arange(P.k, dtype=np.int64)]
    return [np.array(p, dtype=np.int64) for p in itertools.permutations(range(P.k))
            if all(frozenset(p[l] for l in s) == s for s in P.allowed)]


def _compositions(P: Parent, N, LG, conc, include_incomplete):
    k = P.k
    if conc is not None:
        c = np.array(conc, dtype=np.int64)
        return [(c, [p for p in LG if all(c[p[l]] == c[l] for l in range(k))])]
    out = []
    cuts = (itertools.combinations_with_replacement(range(0, N + 1), k - 1) if include_incomplete
            else itertools.combinations(range(1, N), k - 1))
    for cut in cuts:
        c = np.diff((0,) + cut + (N,)).astype(np.int64)
        if not include_incomplete and c.min() == 0:
            continue
        imgs = []
        for p in LG:
            cc = np.empty(k, dtype=np.int64)
            cc[p] = c
            imgs.append(tuple(cc))
        if tuple(c) != max(imgs):
            continue
        out.append((c, [p for p in LG if all(c[p[l]] == c[l] for l in range(k))]))
    return out


def _composition_vector(P: Parent, composition, N):
    if composition is None:
        return None
    if isinstance(composition, dict):
        c = np.zeros(P.k, dtype=np.int64)
        for s, x in composition.items():
            c[P._index(s)] = x
    else:
        c = np.array(composition, dtype=np.int64)
    if len(c) != P.k or c.min() < 0:
        raise ValueError("composition needs one non-negative count per species")
    return c if c.sum() == N else None


# ---------------------------------------------------------------- driver
def _run(P: Parent, n, store, label_exchange, composition, include_incomplete, stats):
    """Enumerate one index n. store=True yields (hnf, labels) per structure; store="batch"
    yields (list of hnfs, labels array) per group of superlattices sharing the same words;
    store=False returns the count."""
    k, d = P.k, P.d
    N = n * d
    conc = _composition_vector(P, composition, N)
    if composition is not None and conc is None:
        return 0
    LG = _label_group(P, label_exchange)
    comps = _compositions(P, N, LG, conc, include_incomplete)
    by_snf = defaultdict(list)
    for H in unique_superlattices(n, P.rots):
        q = Quotient(H)
        by_snf[q.d].append(q)
    total = 0
    for qs in by_snf.values():
        q0 = qs[0]
        phi, A_, m = _column_major(q0, d)
        ident = np.arange(N)

        def conj(p):
            pc = np.empty(N, dtype=np.int64)
            pc[phi] = phi[p]
            return pc
        tp = q0.translation_perms()
        track_perms = []
        for t in range(1, q0.n):
            if q0.elems[t, 2] == 0:
                full = np.concatenate([i * q0.n + tp[t] for i in range(d)])
                track_perms.append(conj(full)[:A_].copy())
        uni, index, members = [], {}, []
        for q in qs:
            mem = set()
            for op in _stabilizer_ops(P, q.H):
                pc = conj(_site_perm(P, q, op))
                if (pc == ident).all():
                    continue
                key = pc.tobytes()
                if key not in index:
                    index[key] = len(uni)
                    uni.append(pc)
                mem.add(index[key])
            members.append(frozenset(mem))
        nu, uni, members = split_reversal(uni, members, m, A_)
        Rp = np.array(uni, dtype=np.int64) if uni else np.zeros((0, N), dtype=np.int64)
        classes = defaultdict(list)
        for q, mem in zip(qs, members):
            classes[mem].append(q)
        keys = list(classes)
        mask = np.zeros((len(keys), max(len(uni), 1)), dtype=np.bool_)
        for ci, key in enumerate(keys):
            mask[ci, list(key)] = True
        q_alpha = k ** A_
        if q_alpha * max(1, len(LG) * (len(track_perms) + 1)) > 4_000_000:
            raise MemoryError(f"letter table too large for Smith normal form {q0.d} (k^{A_} letters)")
        digits = np.array([[(s // k ** (A_ - 1 - t)) % k for t in range(A_)]
                           for s in range(q_alpha)], dtype=np.int64)
        weights = k ** np.arange(A_ - 1, -1, -1)
        symcnt = np.zeros((q_alpha, k), dtype=np.int64)
        for l in range(k):
            symcnt[:, l] = (digits == l).sum(1)
        a = A_ // d
        allowed = np.ones(q_alpha, dtype=np.bool_)
        for t in range(A_):
            allowed &= np.isin(digits[:, t], sorted(P.allowed[t // a]))
        ISP, CM, LT = product_tables(Rp, m, A_, digits, weights)
        for c, lp in comps:
            sw, swt = [], []
            for tpr in [np.arange(A_)] + track_perms:
                for e in lp:
                    is_id_e = (e == np.arange(k)).all()
                    if is_id_e and (tpr == np.arange(A_)).all():
                        continue
                    newd = np.empty_like(digits)
                    newd[:, tpr] = e[digits]
                    sw.append(newd @ weights)
                    swt.append(bool(is_id_e))
            SW = np.array(sw, dtype=np.int64) if sw else np.zeros((0, q_alpha), dtype=np.int64)
            swtrans = np.array(swt, dtype=np.bool_)
            MU = reversal_maps(nu, digits, weights, SW)
            counts = np.zeros(len(keys), dtype=np.int64)
            if not store:
                ns, nodes = torus_kernel(m, A_, k, c, symcnt, SW, swtrans, Rp, mask, counts,
                                         False, np.zeros((0, m), dtype=np.int64), np.zeros(0, dtype=np.int64),
                                         allowed, MU, ISP, CM, LT)
            else:
                # Listing streams: the kernel fills a fixed buffer, pauses, and resumes where it
                # stopped, so memory does not grow with the number of structures.
                cap = max(LIST_BATCH, len(keys))
                ow = np.zeros((cap, m), dtype=np.int64)
                oc = np.zeros(cap, dtype=np.int64)
                W_ = max(SW.shape[0] * m, 1)
                state = (np.zeros(k, dtype=np.int64), np.zeros(m + 1, dtype=np.int64),
                         np.zeros(m + 2, dtype=np.int64), np.full(m + 2, -1, dtype=np.int64),
                         np.zeros((m + 1, W_), dtype=np.int64), np.zeros(m + 1, dtype=np.int64),
                         np.zeros(max(m * MU.shape[0], 1), dtype=np.int64), np.zeros(m + 2, dtype=np.int64),
                         np.zeros(3, dtype=np.int64))
                powers = k ** np.arange(A_ - 1, -1, -1, dtype=np.int64)
                done = False
                while not done:
                    ns, nodes, done = torus_kernel_resume(m, A_, k, c, symcnt, SW, swtrans, Rp, mask, counts,
                                                          ow, oc, allowed, MU, ISP, CM, LT, *state)
                    if ns == 0:
                        continue
                    # decode the batch: letter -> its A_ base-k digits -> site labels
                    flat = ((ow[:ns, :, None] // powers) % k).reshape(ns, m * A_)
                    labs = flat[:, phi]                                # (ns, N), site order
                    cls = oc[:ns].copy()
                    if store == "batch":
                        for ci in np.unique(cls):
                            yield [q.H for q in classes[keys[ci]]], labs[cls == ci]
                    else:
                        for lab_row, ci in zip(labs, cls):
                            lab = tuple(int(x) for x in lab_row)
                            for q in classes[keys[ci]]:
                                yield q.H, lab
            if stats is not None:
                stats["nodes"] = stats.get("nodes", 0) + int(nodes)
            for ci, key in enumerate(keys):
                total += int(counts[ci]) * len(classes[key])
    return total


def _sizes(sizes):
    if isinstance(sizes, (int, np.integer)):
        return [int(sizes)]
    return [int(s) for s in sizes]


def enumerate_structures(parent: Parent, sizes, composition=None, label_exchange=False,
                         include_incomplete=False):
    """Yield every symmetrically distinct derivative structure.

    Parameters
    ----------
    parent : Parent
    sizes : int or iterable of int
        Supercell sizes n (number of primitive cells).
    composition : dict or sequence, optional
        Atom counts per species, e.g. ``{"Cu": 6, "Au": 2}``. Sizes whose site
        count does not match are skipped. Default: all compositions.
    label_exchange : bool
        Also identify structures that differ by a permutation of the species
        (the convention of Hart & Forcade's tables). Default False.
    include_incomplete : bool
        Also list structures in which some species is absent.
    """
    for n in _sizes(sizes):
        yield from (Structure(parent, H, lab) for H, lab in
                    _run(parent, n, True, label_exchange, composition, include_incomplete, None))


def enumerate_labelings(parent: Parent, sizes, composition=None, label_exchange=False,
                        include_incomplete=False):
    """Fast listing: yield the structures as arrays of labelings, without a Python object per structure.

    Yields ``(hnf, labels)`` pairs. ``hnf`` is a superlattice in Hermite normal form and
    ``labels`` is an integer array of shape ``(r, n * d)``: each row is one derivative
    structure on that superlattice, in the site order of :class:`Structure` (``labels[j, i * n + g]``
    is the species index on basis site i at quotient-group element g). Rows come in batches of
    at most ``LIST_BATCH`` words, so memory stays fixed. The same array may be yielded for
    several superlattices that share a Smith normal form and symmetry; treat it as read-only.

    ``Structure(parent, hnf, tuple(row))`` rebuilds the object that :func:`enumerate_structures`
    would give. Parameters are as for :func:`enumerate_structures`.
    """
    dt = np.int8 if parent.k <= 127 else np.int64
    for n in _sizes(sizes):
        for Hs, labs in _run(parent, n, "batch", label_exchange, composition, include_incomplete, None):
            labs = np.ascontiguousarray(labs, dtype=dt)
            labs.flags.writeable = False
            for H in Hs:
                yield H, labs


def count_structures(parent: Parent, sizes, composition=None, label_exchange=False,
                     include_incomplete=False, stats=None):
    """Number of derivative structures for each size, without listing them.

    Returns a dict {n: count}.
    """
    out = {}
    for n in _sizes(sizes):
        gen = _run(parent, n, False, label_exchange, composition, include_incomplete, stats)
        try:
            while True:
                next(gen)
        except StopIteration as stop:
            out[n] = stop.value
    return out
