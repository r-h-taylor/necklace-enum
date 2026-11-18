"""Compiled generator for canonical words (Algorithm 1 of the paper).

A labeling of the quotient group Z_d1 + Z_d2 + Z_m is stored as a word of m
letters, one per column; a letter packs the a entries of its column in base k.
Translations act on words as rotations combined with letterwise maps
(track shifts and label exchanges, table SW). The kernel generates Lyndon words
of fixed content with the FKM recursion, prunes by every letterwise map with
active comparisons, prunes by the reversal-type operation (inversion, tables MU)
with reversed comparisons, and tests the remaining point operations (site
permutations Rp) on each completed word.
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _rot_test(r, word, lab, m, a, k, Rp, SW, img, vsym, ISP, CM, LT):
    """True if some spelling of e.pi.rot(r.w) is smaller than w."""
    n = m * a
    if ISP[r]:
        # product form: column p -> CM[r, p], letter map LT[r]
        for p in range(m):
            vsym[CM[r, p]] = LT[r, word[p]]
    else:
        for f in range(n):
            img[Rp[r, f]] = lab[f]
        for p in range(m):
            s = 0
            for t in range(a):
                s = s * k + img[p * a + t]
            vsym[p] = s
    w0 = word[0]
    nsw = SW.shape[0]
    for i in range(-1, nsw):              # -1 = identity map
        for p in range(m):
            x = vsym[p] if i < 0 else SW[i, vsym[p]]
            if x < w0:
                return True
            if x == w0:
                for d in range(1, m):
                    y = vsym[(p + d) % m]
                    if i >= 0:
                        y = SW[i, y]
                    if y != word[d]:
                        if y < word[d]:
                            return True
                        break
    return False


@njit(cache=True)
def torus_kernel(m, a, k, c, symcnt, SW, swtrans, Rp, cls_mask, counts,
                 store, out_words, out_cls, allowed, MU, ISP, CM, LT):
    # MU: letter maps mu = pi o nu for the reversal-type operation (inversion),
    # shape (nmu, q); nmu = 0 disables reversal pruning.
    q = symcnt.shape[0]
    nmu = MU.shape[0]
    PEND = np.zeros(max(m * nmu, 1), dtype=np.int64)   # pending reversal starts, code = mu*(m+1)+j
    NP = np.zeros(m + 2, dtype=np.int64)
    n = m * a
    nsw = SW.shape[0]
    nr = Rp.shape[0]
    nc = cls_mask.shape[0]
    rem = c.copy()
    A = np.zeros(m + 1, dtype=np.int64)
    P = np.zeros(m + 2, dtype=np.int64)
    S = np.full(m + 2, -1, dtype=np.int64)
    W = max(nsw * m, 1)
    ACT = np.zeros((m + 1, W), dtype=np.int64)   # active (i, j) codes per depth
    NA = np.zeros(m + 1, dtype=np.int64)
    tmp = np.zeros(W, dtype=np.int64)
    word = np.empty(m, dtype=np.int64)
    lab = np.empty(n, dtype=np.int64)
    img = np.empty(n, dtype=np.int64)
    vsym = np.empty(m, dtype=np.int64)
    alive = np.zeros(nc, dtype=np.bool_)
    nstore = 0
    nodes = 0
    t = 1
    P[1] = 1
    while t >= 1:
        if t > m:
            if P[t] == m:                                  # Lyndon under C_m
                ok = True
                for s_ in range(NA[m]):                    # finish wrapped tails
                    code = ACT[m, s_]
                    i = code // m
                    j = code % m
                    res = 0
                    for s in range(1, j + 1):
                        x = SW[i, A[s]]
                        y = A[m - j + s]
                        if x != y:
                            res = -1 if x < y else 1
                            break
                    if res < 0 or (res == 0 and swtrans[i]):
                        ok = False
                        break
                if ok:
                    for s_ in range(NP[m]):
                        code = PEND[s_]
                        mu = code // (m + 1)
                        j = code % (m + 1)
                        for ii in range(j + 1, m + 1):
                            x = MU[mu, A[m - (ii - j - 1)]]
                            y = A[ii]
                            if x != y:
                                if x < y:
                                    ok = False
                                break
                        if not ok:
                            break
                if ok:
                    for p in range(m):
                        word[p] = A[p + 1]
                        s = A[p + 1]
                        for tt in range(a - 1, -1, -1):
                            lab[p * a + tt] = s % k
                            s //= k
                    for cl in range(nc):
                        alive[cl] = True
                    nalive = nc
                    for r in range(nr):
                        need = False
                        for cl in range(nc):
                            if alive[cl] and cls_mask[cl, r]:
                                need = True
                                break
                        if not need:
                            continue
                        if _rot_test(r, word, lab, m, a, k, Rp, SW, img, vsym, ISP, CM, LT):
                            for cl in range(nc):
                                if cls_mask[cl, r] and alive[cl]:
                                    alive[cl] = False
                                    nalive -= 1
                            if nalive == 0:
                                break
                    for cl in range(nc):
                        if alive[cl]:
                            counts[cl] += 1
                            if store and nstore < out_words.shape[0]:
                                for p in range(m):
                                    out_words[nstore, p] = word[p]
                                out_cls[nstore] = cl
                                nstore += 1
            t -= 1
            continue
        nodes += 1
        lo = A[t - P[t]]
        if S[t] == -1:
            cand = lo
        else:
            for l in range(k):
                rem[l] += symcnt[S[t], l]
            cand = S[t] + 1
        chosen = -1
        na = 0
        while cand < q:
            ok = allowed[cand]           # site restrictions: letter must be allowed
            for l in range(k):
                if not ok:
                    break
                if symcnt[cand, l] > rem[l]:
                    ok = False
                    break
            if ok:
                na = 0
                # carried-over comparisons, then the ones starting at j = t-1
                for s_ in range(NA[t - 1] + nsw):
                    if s_ < NA[t - 1]:
                        code = ACT[t - 1, s_]
                    else:
                        code = (s_ - NA[t - 1]) * m + (t - 1)
                    i = code // m
                    j = code % m
                    x = SW[i, cand]
                    y = cand if j == 0 else A[t - j]
                    if x < y:
                        ok = False
                        break
                    if x == y:
                        tmp[na] = code
                        na += 1
            if ok and nmu > 0:
                npend = NP[t - 1]
                for mu in range(nmu):
                    tie = True
                    for ii in range(1, t + 1):
                        src = t - ii + 1
                        xl = cand if src == t else A[src]
                        x = MU[mu, xl]
                        y = cand if ii == t else A[ii]
                        if x != y:
                            tie = False
                            if x < y:
                                ok = False
                            break
                    if not ok:
                        break
                    if tie:
                        PEND[npend] = mu * (m + 1) + t
                        npend += 1
                if ok:
                    NP[t] = npend
            elif ok:
                NP[t] = NP[t - 1]
            if ok:
                chosen = cand
                break
            cand += 1
        if chosen == -1:
            S[t] = -1
            t -= 1
            continue
        S[t] = chosen
        A[t] = chosen
        for l in range(k):
            rem[l] -= symcnt[chosen, l]
        for s_ in range(na):
            ACT[t, s_] = tmp[s_]
        NA[t] = na
        P[t + 1] = P[t] if chosen == lo else t
        t += 1
        if t <= m:
            S[t] = -1
    return nstore, nodes


def split_reversal(uni, members, m, a):
    """Find an operation acting on words as reversal followed by a letterwise
    track permutation nu (p -> -p, track t -> nu[t]) that every superlattice of
    the SNF has. Returns (nu or None, uni', members') with it removed."""
    if not uni:
        return None, uni, members
    common = frozenset.intersection(*members) if members else frozenset()
    cols = ((-np.arange(m)) % m) * a
    for r in sorted(common):
        pc = uni[r]
        nu = pc[:a]
        if (nu >= a).any():
            continue
        if (pc.reshape(m, a) == cols[:, None] + nu[None, :]).all():
            keep = [i for i in range(len(uni)) if i != r]
            remap = {old: new for new, old in enumerate(keep)}
            uni2 = [uni[i] for i in keep]
            mem2 = [frozenset(remap[i] for i in mem if i != r) for mem in members]
            return nu.copy(), uni2, mem2
    return None, uni, members


def product_tables(Rp, m, a, digits, weights):
    """For each rotation permutation of the form (p,t) -> (CM[p], tm[t]),
    return ISP=True, the column map CM and the letter table LT."""
    nr = Rp.shape[0]
    q_alpha = digits.shape[0]
    ISP = np.zeros(nr, dtype=np.bool_)
    CM = np.zeros((max(nr, 1), m), dtype=np.int64)
    LT = np.zeros((max(nr, 1), q_alpha), dtype=np.int64)
    for r in range(nr):
        P = Rp[r].reshape(m, a)
        cols = P // a
        tr = P % a
        if (cols == cols[:, :1]).all() and (tr == tr[:1]).all():
            ISP[r] = True
            CM[r] = cols[:, 0]
            newd = np.empty_like(digits)
            newd[:, tr[0]] = digits
            LT[r] = newd @ weights
    return ISP, CM, LT


def reversal_maps(nu, digits, weights, SW):
    """Letter tables for mu = nu and mu = pi o nu over all letterwise maps pi."""
    q_alpha = digits.shape[0]
    if nu is None:
        return np.zeros((0, q_alpha), dtype=np.int64)
    newd = np.empty_like(digits)
    newd[:, nu] = digits
    NU = newd @ weights
    rows = [NU] + [SW[i][NU] for i in range(SW.shape[0])]
    return np.array(rows, dtype=np.int64)
