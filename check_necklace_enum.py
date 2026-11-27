#!/usr/bin/env python3
"""Stand-alone check of necklace-enum.

Run from the package folder:

    python check_necklace_enum.py          # under a minute
    python check_necklace_enum.py --full   # more cases and larger cells (a few minutes)

Section 4 checks the counts against our own reimplementation of the table-based Hart & Forcade
(2008) algorithm (tests/hart_forcade.py). The reimplementation is not enumlib.

Uses the installed package if there is one, otherwise the copy in ./src.
Needs numpy and numba; ase and spglib are optional (export check).
"""
import argparse
import itertools
import os
import sys
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import necklace_enum as te
except ImportError:
    sys.path.insert(0, os.path.join(HERE, "src"))
    import necklace_enum as te
sys.path.insert(0, os.path.join(HERE, "tests"))
from necklace_enum.lattice import unique_superlattices  # noqa: E402

RESULTS = []


def report(name, ok, detail=""):
    RESULTS.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


# ---------------------------------------------------------------- 1. published counts
# Hart & Forcade, Phys. Rev. B 77, 224115 (2008): Tables IV-VII (label exchange, every species present)
FCC_SUPER = [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]
FCC2 = [2, 3, 12, 14, 50, 52, 229, 252, 685, 682, 3875, 2624, 9628, 16584, 49764, 42135,
        212612, 174104, 867893, 1120708, 2628180, 3042732]
FCC3 = {3: 3, 4: 13, 5: 23, 6: 130, 7: 197, 8: 1267, 9: 2322, 10: 9332}
FCC4 = {4: 7, 5: 9, 6: 110, 7: 211, 8: 2110, 9: 5471, 10: 32362}
SC2 = {2: 3, 3: 3, 4: 15}


def check_published(full):
    print("1. Counts published by Hart & Forcade (2008)")
    got = [len(unique_superlattices(n, te.fcc().rots)) for n in range(1, 11)]
    report("fcc superlattices, n = 1-10 (Table IV)", got == FCC_SUPER)
    nmax = 23 if full else 17
    c = te.count_structures(te.fcc(2), range(2, nmax + 1), label_exchange=True)
    report(f"binary fcc, n = 2-{nmax} (Table V)", [c[n] for n in range(2, nmax + 1)] == FCC2[:nmax - 1])
    report("ternary fcc, n = 3-10 (Table VII)",
           te.count_structures(te.fcc(3), list(FCC3), label_exchange=True) == FCC3)
    report("quaternary fcc, n = 4-10 (Table VII)",
           te.count_structures(te.fcc(4), list(FCC4), label_exchange=True) == FCC4)
    report("binary simple cubic, n = 2-4 (Table VI)",
           te.count_structures(te.sc(2), list(SC2), label_exchange=True) == SC2)


# ---------------------------------------------------------------- 2. brute force
def check_brute_force(full):
    print("2. Structure-by-structure comparison with an independent brute-force code")
    from brute import brute, canon, check_superlattices, label_group, to_cell_labeling
    parents = {"fcc binary": te.fcc(2), "fcc ternary": te.fcc(3), "bcc binary": te.bcc(2),
               "sc ternary": te.sc(3), "hcp binary": te.hcp(2), "hcp ternary": te.hcp(3),
               "hcp restricted {A,B}|{B,C}": te.hcp(3, allowed=[{0, 1}, {1, 2}])}
    nmax = {"fcc binary": 7 if full else 5, "fcc ternary": 4 if full else 3,
            "bcc binary": 6 if full else 4, "sc ternary": 3, "hcp binary": 4 if full else 3,
            "hcp ternary": 2, "hcp restricted {A,B}|{B,C}": 3}
    cases = [(name, n, le, None) for name in parents for n in range(1, nmax[name] + 1) for le in (False, True)]
    cases += [("fcc binary", 8, False, (5, 3)), ("hcp binary", 4, False, (5, 3)),
              ("hcp restricted {A,B}|{B,C}", 4, True, (3, 2, 3))]
    for name, n, le, conc in cases:
        P = parents[name]
        reps = unique_superlattices(n, P.rots)
        check_superlattices(P, n, reps)
        ref = brute(P, n, reps, le, conc)
        LG = label_group(P, le)
        got = defaultdict(list)
        for s in te.enumerate_structures(P, n, composition=conc, label_exchange=le):
            cell, group, _ = ref[tuple(s.hnf.ravel())]
            got[tuple(s.hnf.ravel())].append(canon(to_cell_labeling(P, cell, s.hnf, s.labels), group, LG))
        ok = all(len(got.get(h, [])) == len(set(got.get(h, []))) and set(got.get(h, [])) == e
                 for h, (_, _, e) in ref.items())
        total = sum(len(e) for _, _, e in ref.values())
        ok = ok and te.count_structures(P, n, composition=conc, label_exchange=le)[n] == total
        tag = f"{name}, n={n}" + (", label exchange" if le else "") + (f", composition {conc}" if conc else "")
        report(tag, ok, f"{total} structures")


# ---------------------------------------------------------------- 3. export
def check_export(full):
    print("3. Exported crystals are primitive in their own cells and pairwise distinct")
    try:
        import spglib
        from ase.utils.structure_comparator import SymmetryEquivalenceCheck
    except ImportError:
        print("  [SKIP] needs ase and spglib (pip install ase spglib)")
        return
    comp = SymmetryEquivalenceCheck()
    for P, sizes, label in [(te.fcc(["Cu", "Au"]), range(1, 8 if full else 6), "fcc Cu-Au"),
                            (te.hcp(["Mg", "Cd"]), range(1, 4), "hcp Mg-Cd")]:
        structs = list(te.enumerate_structures(P, sizes))
        atoms = [s.to_ase() for s in structs]
        prim = all(len(spglib.find_primitive((a.cell[:], a.get_scaled_positions(), a.get_atomic_numbers()),
                                             symprec=1e-5)[2]) == len(a) for a in atoms)
        # candidates for equivalence must share formula and space group (ASE's comparator alone
        # can report false matches between structures of different symmetry)
        by_formula = defaultdict(list)
        for a in atoms:
            sg = spglib.get_symmetry_dataset((a.cell[:], a.get_scaled_positions(), a.get_atomic_numbers()),
                                             symprec=1e-5).number
            by_formula[(a.get_chemical_formula(), sg)].append(a)
        distinct = all(not comp.compare(x, y) for g in by_formula.values()
                       for x, y in itertools.combinations(g, 2))
        report(f"{label}, sizes {sizes.start}-{sizes.stop - 1}", prim and distinct,
               f"{len(structs)} structures")


# ---------------------------------------------------------------- 4. counts against the 2008 algorithm (reimplementation)
def check_reference(full):
    print("4. Counts against our reimplementation of the table-based Hart & Forcade (2008) algorithm")
    print("   (tests/hart_forcade.py; not enumlib; all compositions, label exchange)")
    import hart_forcade as hf
    cases = [(2, n) for n in (16, 18, 20, 22, 24)] + [(3, n) for n in (10, 12, 14)]
    if full:
        cases = [(2, n) for n in (16, 18, 20, 22, 24, 25, 26, 28)] + [(3, n) for n in (10, 12, 14, 16)]
    for k, n in cases:
        P = te.fcc(k)
        name = f"fcc {'binary' if k == 2 else 'ternary'}, n={n}"
        if 8 * k ** n > 2e9:
            continue                                      # the reference needs 8 k^n bytes
        ours = te.count_structures(P, n, label_exchange=True)[n]
        ref = hf.count_structures(P, n)
        report(name, ref == ours, f"{ours:,d} structures" + ("" if ref == ours else f", reference {ref:,d}"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", action="store_true", help="more cases and larger cells")
    args = ap.parse_args()
    print(f"necklace-enum {te.__version__} from {os.path.dirname(te.__file__)}\n")
    t0 = time.perf_counter()
    for step in (check_published, check_brute_force, check_export, check_reference):
        step(args.full)
        print()
    failed = [n for n, ok in RESULTS if not ok]
    print(f"{len(RESULTS) - len(failed)} of {len(RESULTS)} checks passed in {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED:", *failed, sep="\n  ")
    sys.exit(1 if failed else 0)
