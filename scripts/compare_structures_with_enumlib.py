"""Structure-by-structure comparison of necklace-enum with enumlib.

    python3 scripts/compare_structures_with_enumlib.py ~/Projects/enumlib/src/enum.x
    python3 scripts/compare_structures_with_enumlib.py ~/Projects/enumlib/src/enum.x --makestr ~/Projects/enumlib/aux_src/makeStr.py

For each case both programs run on the same struct_enum.in. Both outputs are turned into
crystals by enumlib's own makeStr.py (so necklace-enum's struct_enum.out is read exactly as
enumlib's is), and every crystal is reduced to a canonical key by code that shares nothing with
either enumerator: the cell is brought to a fixed representative of its superlattice class by a
parent symmetry operation, and the labeling is minimised over the full group of that superlattice
(translations x the parent operations that fix it x label permutations), using the independent
brute-force module in tests/brute.py. The two sets of keys must be equal, with no duplicates.

Conventions are enumlib's: a "part" listing identifies structures related by a permutation of the
species and has every species present; a "full" listing does neither; fixed concentrations
are listed without label exchange.
"""
import argparse
import collections
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tests"))
import necklace_enum as te                                       # noqa: E402
from necklace_enum.lattice import unique_superlattices          # noqa: E402
from brute import Cell, canon, check_superlattices, is_int, label_group   # noqa: E402

# name, parent factory, species k, sizes, mode ("part" / "full"), concentration (counts per cell) or None
CASES = [
    ("fcc binary", te.fcc, 2, range(1, 11), "part", None),
    ("fcc binary", te.fcc, 2, range(1, 9), "full", None),
    ("fcc ternary", te.fcc, 3, range(1, 7), "part", None),
    ("fcc ternary", te.fcc, 3, range(1, 6), "full", None),
    ("fcc quaternary", te.fcc, 4, range(1, 6), "part", None),
    ("bcc binary", te.bcc, 2, range(1, 10), "part", None),
    ("sc binary", te.sc, 2, range(1, 9), "part", None),
    ("hcp binary", te.hcp, 2, range(1, 6), "part", None),
    ("hcp binary", te.hcp, 2, range(1, 5), "full", None),
    ("hcp ternary", te.hcp, 3, range(1, 4), "part", None),
    ("fcc 9:3", te.fcc, 2, [12], "full", (9, 3)),
    ("fcc 6:3:3", te.fcc, 3, [12], "full", (6, 3, 3)),
    ("bcc 7:3", te.bcc, 2, [10], "full", (7, 3)),
    ("hcp 8:4", te.hcp, 2, [6], "full", (8, 4)),
]
QUICK = [CASES[0][:3] + (range(1, 7),) + CASES[0][4:], CASES[7][:3] + (range(1, 4),) + CASES[7][4:],
         ("fcc 6:2", te.fcc, 2, [8], "full", (6, 2))]


def struct_enum_in(P, sizes, mode, conc):
    k = P.k
    lab = "/".join(map(str, range(k)))
    L = [f"necklacecmp\nbulk"]
    for c in range(3):
        L.append(" ".join(f"{x:.12f}" for x in P.A[:, c]))
    L.append(f" {k} -nary case")
    L.append(f" {P.d}")
    for b in P.basis:
        L.append(" ".join(f"{x:.12f}" for x in b) + f" {lab}")
    L.append(f" {min(sizes)} {max(sizes)}")
    L.append("0.10000000E-06")
    L.append(f"{mode} list of labelings")
    if conc:
        N = max(sizes) * P.d
        L.append("# Concentration restrictions")
        L += [f"{c} {c} {N}" for c in conc]
    return "\n".join(L) + "\n"


def read_poscar(path):
    lines = open(path).read().split("\n")
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[2 + i].split()[:3]] for i in range(3)]) * scale
    counts = [int(x) for x in lines[5].split()]
    assert lines[6].strip().lower().startswith("d"), path
    frac = np.array([[float(x) for x in lines[7 + i].split()[:3]] for i in range(sum(counts))])
    labels = [l for l, c in enumerate(counts) for _ in range(c)]
    return lat, frac @ lat, labels


class Keyer:
    """Canonical key of a crystal on parent P, independent of both enumerators."""

    def __init__(self, P, le):
        self.P, self.LG, self.reps, self.cells = P, label_group(P, le), {}, {}

    def reps_for(self, n):
        if n not in self.reps:
            reps = unique_superlattices(n, self.P.rots)
            check_superlattices(self.P, n, reps)          # independent pairwise test
            self.reps[n] = [(H, np.linalg.inv(H)) for H in reps]
        return self.reps[n]

    def key(self, lat, cart, labels):
        P = self.P
        S = P.Ainv @ lat.T
        assert is_int(S), "cell is not a superlattice of the parent"
        S = np.round(S).astype(int)
        n = round(abs(np.linalg.det(S)))
        for M, tau, _, _ in P.ops:
            S2 = M @ S
            for H, Hi in self.reps_for(n):
                if is_int(Hi @ S2) and round(abs(np.linalg.det(H))) == n:
                    hk = tuple(H.ravel())
                    if hk not in self.cells:
                        cell = Cell(P, H)
                        self.cells[hk] = (cell, cell.group()[0])
                    cell, group = self.cells[hk]
                    R = P.A @ M @ P.Ainv
                    lab = [None] * (n * P.d)
                    for x, l in zip(cart, labels):
                        y = R @ x + tau
                        for j in range(P.d):
                            f = P.Ainv @ (y - P.basis[j])
                            if is_int(f):
                                lab[cell.site(j, np.round(f).astype(int))] = l
                                break
                        else:
                            raise AssertionError("atom not on a parent site")
                    assert None not in lab, "two atoms on one site"
                    return (n, hk, canon(lab, group, self.LG))
        raise AssertionError("no superlattice representative found")


def run_case(enumx, makestr, name, factory, k, sizes, mode, conc, keep):
    P = factory(k)
    le = mode == "part" and conc is None
    work = tempfile.mkdtemp(prefix="cmp_")
    text = struct_enum_in(P, list(sizes), mode, conc)
    keys = {}
    times = {}
    for prog in ("enumlib", "necklace-enum"):
        d = os.path.join(work, prog)
        os.makedirs(d)
        open(os.path.join(d, "struct_enum.in"), "w").write(text)
        t = time.perf_counter()
        if prog == "enumlib":
            subprocess.run([enumx], cwd=d, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            te.run_enumlib_input(os.path.join(d, "struct_enum.in"), os.path.join(d, "struct_enum.out"))
        times[prog] = time.perf_counter() - t
        p = subprocess.run([sys.executable, makestr, "all", "-input", "struct_enum.out", "-outfile", "vasp.{}"],
                           cwd=d, capture_output=True, text=True)
        if p.returncode:
            raise RuntimeError(f"makeStr.py failed for {prog}: {p.stderr.strip()[-300:]}")
        kr = Keyer(P, le)
        keys[prog] = collections.Counter(kr.key(*read_poscar(f)) for f in glob.glob(os.path.join(d, "vasp.*")))
    a, b = keys["enumlib"], keys["necklace-enum"]
    dup_a = sum(c - 1 for c in a.values() if c > 1)
    dup_b = sum(c - 1 for c in b.values() if c > 1)
    only_a, only_b = set(a) - set(b), set(b) - set(a)
    ok = not (dup_a or dup_b or only_a or only_b)
    label = f"{name}, {mode}" + (f", sizes {min(sizes)}-{max(sizes)}" if len(sizes) > 1 else f", n={sizes[0]}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:42s} enumlib {sum(a.values()):6d}  necklace-enum {sum(b.values()):6d}"
          + ("" if ok else f"  | only enumlib {len(only_a)}, only necklace {len(only_b)}, duplicates {dup_a}/{dup_b}"),
          flush=True)
    if not ok:
        print(f"         kept for inspection: {work}")
    elif not keep:
        shutil.rmtree(work)
    return ok, sum(b.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("enumx")
    ap.add_argument("--makestr", default=None, help="enumlib's aux_src/makeStr.py (default: next to enum.x)")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--keep", action="store_true", help="keep the working folders")
    a = ap.parse_args()
    enumx = os.path.abspath(os.path.expanduser(a.enumx))
    makestr = os.path.abspath(os.path.expanduser(a.makestr or os.path.join(os.path.dirname(enumx), "..", "aux_src", "makeStr.py")))
    for f in (enumx, makestr):
        if not os.path.exists(f):
            sys.exit(f"not found: {f}")
    print(f"enumlib: {enumx}\nmakeStr: {makestr}\n")
    t = time.perf_counter()
    res = [run_case(enumx, makestr, *c, keep=a.keep) for c in (QUICK if a.quick else CASES)]
    n_ok = sum(ok for ok, _ in res)
    print(f"\n{n_ok} of {len(res)} cases identical, {sum(n for _, n in res):,d} structures compared, "
          f"{time.perf_counter() - t:.0f} s")
    sys.exit(0 if n_ok == len(res) else 1)


if __name__ == "__main__":
    main()
