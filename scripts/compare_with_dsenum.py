"""Compare necklace-enum with dsenum (Shinohara et al. 2020, ZDD method) on identical problems.

    pip install dsenum
    python scripts/compare_with_dsenum.py            # default problem set
    python scripts/compare_with_dsenum.py --quick    # small sizes, for checking the setup

dsenum does not remove label exchange, so both codes run without it. Both exclude
superperiodic structures and structures that do not use every species. For each problem
it reports the count from each code (they must agree) and three times:
  count: dsenum ZddStructureEnumerator.count()  vs necklace_enum.count_structures
  list:  dsenum, building every labeling from the ZDD for every superlattice
         vs necklace_enum.enumerate_structures, consuming every structure
Each time is the best of --repeat runs, after one warm-up run (Numba compilation is
excluded). Times include all per-superlattice setup in both codes.
"""
import argparse
import platform
import time

import numpy as np

import necklace_enum as te
from dsenum import ZddStructureEnumerator
from dsenum.permutation_group import DerivativeStructurePermutation
from pymatgen.core import Lattice, Structure


def fcc_base():
    return Structure(Lattice([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]), ["Cu"], [[0, 0, 0]])


def hcp_base():
    a, c = 1.0, np.sqrt(8 / 3)
    lat = Lattice.hexagonal(a, c)
    return Structure(lat, ["Mg", "Mg"], [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]])


# name, necklace-enum parent factory, dsenum base, index, k, composition (counts per cell) or None
PROBLEMS = [
    ("fcc binary, all compositions, n=16", te.fcc, fcc_base, 16, 2, None),
    ("fcc binary, all compositions, n=20", te.fcc, fcc_base, 20, 2, None),
    ("fcc binary, all compositions, n=24", te.fcc, fcc_base, 24, 2, None),
    ("fcc ternary, all compositions, n=12", te.fcc, fcc_base, 12, 3, None),
    ("fcc ternary, all compositions, n=14", te.fcc, fcc_base, 14, 3, None),
    ("fcc quaternary, all compositions, n=10", te.fcc, fcc_base, 10, 4, None),
    ("fcc 15:5, n=20", te.fcc, fcc_base, 20, 2, (15, 5)),
    ("fcc 24:8, n=32", te.fcc, fcc_base, 32, 2, (24, 8)),
    ("fcc 10:5:5, n=20", te.fcc, fcc_base, 20, 3, (10, 5, 5)),
    ("hcp binary, all compositions, 24 sites", te.hcp, hcp_base, 12, 2, None),
    ("hcp ternary, all compositions, 16 sites", te.hcp, hcp_base, 8, 3, None),
]
QUICK = [
    ("fcc binary, all compositions, n=8", te.fcc, fcc_base, 8, 2, None),
    ("fcc ternary, all compositions, n=6", te.fcc, fcc_base, 6, 3, None),
    ("fcc 6:2, n=8", te.fcc, fcc_base, 8, 2, (6, 2)),
    ("hcp binary, all compositions, 12 sites", te.hcp, hcp_base, 6, 2, None),
]


def best(f, repeat):
    f()                                                    # warm-up / JIT
    ts, out = [], None
    for _ in range(repeat):
        t = time.perf_counter()
        out = f()
        ts.append(time.perf_counter() - t)
    return out, min(ts)


def ds_enum(base, n, k, comp):
    ratio = None
    if comp is not None:
        g = np.gcd.reduce(comp)
        ratio = [c // g for c in comp]
    return ZddStructureEnumerator(base, n, k, composition_constraints=ratio,
                                  remove_superperiodic=True, remove_incomplete=True, verbose=False)


def ds_count(base, n, k, comp):
    return ds_enum(base, n, k, comp).count()


def ds_list(base, n, k, comp):
    e = ds_enum(base, n, k, comp)
    disp = e.base_structure.frac_coords
    total = 0
    for hnf in e.list_reduced_HNF:
        perm = DerivativeStructurePermutation(hnf, disp, e.rotations, e.translations)
        for _ in e._generate_coloring_with_hnf(hnf, perm, [], []):
            total += 1
    return total


def ne_count(P, n, comp):
    return te.count_structures(P, n, composition=comp, label_exchange=False)[n]


def ne_list(P, n, comp):
    return sum(1 for _ in te.enumerate_structures(P, n, composition=comp, label_exchange=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--no-list", action="store_true", help="skip the listing comparison")
    a = ap.parse_args()
    print(f"# {platform.platform()}, {platform.processor() or platform.machine()}, Python {platform.python_version()}")
    print(f"{'problem':42s} {'structures':>12s} {'count: dsenum':>14s} {'necklace':>9s}"
          + ("" if a.no_list else f" {'list: dsenum':>13s} {'necklace':>9s}"))
    for name, parent, base, n, k, comp in (QUICK if a.quick else PROBLEMS):
        P, B = parent(k), base()
        c_ds, t_ds = best(lambda: ds_count(B, n, k, comp), a.repeat)
        c_ne, t_ne = best(lambda: ne_count(P, n, comp), a.repeat)
        line = f"{name:42s} {c_ne:>12,d} {t_ds:14.2f} {t_ne:9.2f}"
        if c_ds != c_ne:
            line += f"  COUNT MISMATCH: dsenum {c_ds:,d}"
        if not a.no_list:
            l_ds, tl_ds = best(lambda: ds_list(B, n, k, comp), 1)
            l_ne, tl_ne = best(lambda: ne_list(P, n, comp), 1)
            line += f" {tl_ds:13.2f} {tl_ne:9.2f}"
            if not (l_ds == l_ne == c_ne):
                line += f"  LIST MISMATCH: {l_ds:,d} / {l_ne:,d}"
        print(line, flush=True)


if __name__ == "__main__":
    main()
