"""struct_enum.out output follows enumlib's conventions.

These tests do not need enumlib. Crystals are rebuilt from the written file with
enumlib's indexing rule (as in its makeStr.py: site at lattice point z has address
g = L z mod D and label index i*n + g1*d2*d3 + g2*d3 + g3) and compared with the
structures necklace-enum produces directly. Degeneracy columns are checked
against a direct count of (superlattice, labeling) pairs.
"""
import io
import itertools

import numpy as np
import pytest

import necklace_enum as te
from necklace_enum.lattice import Quotient, hnfs

ENUMLIB_EXAMPLE_IN = """fcc
bulk
0.50000000     0.50000000      0.0000000        # a1 parent lattice vector
0.50000000      0.0000000     0.50000000        # a2 parent lattice vector
 0.0000000     0.50000000     0.50000000        # a3 parent lattice vector
 2 -nary case
    1 # Number of points in the multilattice
 0.0000000      0.0000000      0.0000000    0/1   # d01 d-vector
    1 8   # Starting and ending cell sizes for search
0.10000000E-06 # Epsilon (finite precision parameter)
full list of labelings
# Concentration restrictions
1 2 2  # Or you could use 50 100 100
0 2 2  # Or you could use 0 100 100
"""


def rebuild(header, row):
    """Labels at every lattice point of the supercell, keyed by (site kind, fractional
    position in the supercell), using enumlib's convention."""
    A = header["lattice"].T
    H, L, D = row["hnf"], row["L"], row["snf"]
    n = row["n"]
    Hi = np.linalg.inv(H)
    out = {}
    nd = len(header["basis"])
    r = int(np.abs(H).max()) + 1
    for z in itertools.product(range(-r, r + 1), repeat=3):
        z = np.array(z)
        f = Hi @ z
        f = tuple(np.round(f - np.floor(f + 1e-9), 6) % 1.0)
        if (0, f) in out:
            continue
        g = (L @ z) % D
        gi = g[0] * D[1] * D[2] + g[1] * D[2] + g[2]
        for i in range(nd):
            out[(i, f)] = int(row["labeling"][i * n + gi])
    return out


def from_structure(s):
    q = Quotient(s.hnf)
    Hi = np.linalg.inv(s.hnf)
    out = {}
    for i in range(s.parent.d):
        for g in range(q.n):
            z = q.Linv @ q.elems[g]
            f = Hi @ z
            out[(i, tuple(np.round(f - np.floor(f + 1e-9), 6) % 1.0))] = s.labels[i * q.n + g]
    return out


@pytest.mark.parametrize("parent,sizes,kw", [
    (te.fcc(2), range(1, 7), dict(include_incomplete=True)),
    (te.hcp(2), range(1, 4), dict(include_incomplete=True)),
    (te.hcp(3, allowed=[{0, 1}, {1, 2}]), range(1, 4), dict(include_incomplete=True)),
    (te.fcc(2), 8, dict(composition=(6, 2))),
])
def test_rows_rebuild_the_same_structures(parent, sizes, kw):
    buf = io.StringIO()
    n_written = te.write_struct_enum(buf, parent, sizes, **kw)
    buf.seek(0)
    header, rows = te.read_struct_enum(buf)
    structs = list(te.enumerate_structures(parent, sizes, **kw))
    assert n_written == len(rows) == len(structs)
    by_key = {(tuple(s.hnf.ravel()), "".join(map(str, s.labels))): s for s in structs}
    assert len(by_key) == len(structs)
    for row in rows:
        s = by_key[(tuple(row["hnf"].ravel()), row["labeling"])]
        assert rebuild(header, row) == from_structure(s)


def _direct_pairs(P, n, include_incomplete):
    total = 0
    for H in hnfs(n):
        q = Quotient(H)
        tp = q.translation_perms()
        shifts = [np.concatenate([i * n + tp[t] for i in range(P.d)]) for t in range(1, n)]
        for lab in itertools.product(*[sorted(P.allowed[s // n]) for s in range(n * P.d)]):
            lab = np.array(lab)
            if not include_incomplete and np.bincount(lab, minlength=P.k).min() == 0:
                continue
            if any((lab[p] == lab).all() for p in shifts):
                continue
            total += 1
    return total


@pytest.mark.parametrize("parent,nmax,include_incomplete", [
    (te.fcc(2), 6, True), (te.fcc(2), 6, False), (te.fcc(3), 4, False), (te.hcp(2), 3, True)])
def test_degeneracies_count_every_labeling(parent, nmax, include_incomplete):
    buf = io.StringIO()
    te.write_struct_enum(buf, parent, range(1, nmax + 1), include_incomplete=include_incomplete)
    buf.seek(0)
    _, rows = te.read_struct_enum(buf)
    for n in range(1, nmax + 1):
        mine = sum(r["tot_degen"] for r in rows if r["n"] == n)
        assert mine == _direct_pairs(parent, n, include_incomplete)
        assert all(r["tot_degen"] == r["hnf_degen"] * r["lab_degen"] for r in rows)


def test_partial_listing_uses_label_exchange_like_enumlib():
    # enumlib's "partial" listing identifies structures related by swapping species,
    # so its counts are Hart & Forcade's label-exchange counts.
    buf = io.StringIO()
    te.write_struct_enum(buf, te.fcc(2), range(1, 9))
    buf.seek(0)
    _, rows = te.read_struct_enum(buf)
    counts = [sum(1 for r in rows if r["n"] == n) for n in range(2, 9)]
    assert counts == [2, 3, 12, 14, 50, 52, 229]


def test_reads_enumlib_example_input(tmp_path):
    spec = te.read_struct_enum_in(io.StringIO(ENUMLIB_EXAMPLE_IN))
    assert spec["sizes"] == list(range(1, 9))
    assert spec["concentration_ranges"] == [(1, 2, 2), (0, 2, 2)]
    assert spec["include_incomplete"] is True
    (tmp_path / "struct_enum.in").write_text(ENUMLIB_EXAMPLE_IN)
    n = te.run_enumlib_input(tmp_path / "struct_enum.in", tmp_path / "struct_enum.out")
    assert n == 376                                     # enumlib's enum.x writes 376 for this input
    _, rows = te.read_struct_enum(tmp_path / "struct_enum.out")
    assert all(r["labeling"].count("0") * 2 >= len(r["labeling"]) for r in rows)
