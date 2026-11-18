"""Reproduce the numbers quoted in the paper (this code's columns only).

    python scripts/reproduce_paper.py          # quick: counts, cyclic statistics, Table 1
    python scripts/reproduce_paper.py --full   # also binary fcc n = 27..32 (about ten minutes)

enumlib timings in Table 1 require a separate enumlib build; see the paper.
"""
import argparse
import time

import necklace_enum as te
from necklace_enum.lattice import Quotient, unique_superlattices


def cyclic_stats():
    P = te.fcc(2)
    tot = cyc = 0
    for n in range(2, 33):
        for H in unique_superlattices(n, P.rots):
            tot += 1
            cyc += Quotient(H).d[:2] == (1, 1)
    print(f"fcc superlattices n=2..32: {tot}, cyclic quotient group: {cyc}")
    reps = unique_superlattices(24, P.rots)
    cyc24 = [H for H in reps if Quotient(H).d[:2] == (1, 1)]
    print(f"n=24: {len(reps)} superlattices, {len(cyc24)} cyclic")


def table1():
    rows = [("fcc binary, all compositions, n=24", te.fcc(2), 24, None, True),
            ("fcc binary, all compositions, n=26", te.fcc(2), 26, None, True),
            ("fcc ternary, all compositions, n=16", te.fcc(3), 16, None, True),
            ("fcc quaternary, all compositions, n=12", te.fcc(4), 12, None, True),
            ("fcc 10:5:5, n=20", te.fcc(3), 20, (10, 5, 5), False),
            ("fcc 24:8, n=32", te.fcc(2), 32, (24, 8), False),
            ("hcp ternary, all compositions, 16 sites", te.hcp(3), 8, None, True),
            ("hcp 12:12, 24 sites", te.hcp(2), 12, (12, 12), False)]
    te.count_structures(te.fcc(2), 4)                       # compile
    print(f"{'problem':42s} {'structures':>14s} {'time (s)':>9s}")
    for name, P, n, comp, le in rows:
        t = time.perf_counter()
        c = te.count_structures(P, n, composition=comp, label_exchange=le)[n]
        print(f"{name:42s} {c:>14,d} {time.perf_counter() - t:9.2f}")


def binary_fcc(nmax):
    for n in range(2, nmax + 1):
        t = time.perf_counter()
        c = te.count_structures(te.fcc(2), n, label_exchange=True)[n]
        print(f"binary fcc n={n:2d}: {c:>14,d}  ({time.perf_counter() - t:.1f} s)", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    cyclic_stats()
    table1()
    binary_fcc(32 if a.full else 26)
