"""Time necklace-enum alone on the problems of the paper's Table 1: counting, listing as Python
objects, and fast listing as arrays of labelings.

    python3 scripts/time_necklace.py --cooldown 60 | tee scripts/necklace_timings_M2.txt

All runs are without label exchange and exclude structures that lack an atom type (the convention
of Table 1). Each problem gets a warm-up run first (Numba compiles on first use). Counting is the
best of three runs; each listing is one run. --cooldown pauses before every timed run, so that a
fanless machine starts each run cool.

  count   count_structures: the generator alone, nothing handed to Python
  arrays  enumerate_labelings(copy=True): every structure handed to Python as its own row of an
          integer array
  objects enumerate_structures(degeneracy=False): every structure handed to Python as a Structure
          object, without degeneracies
"""
import argparse
import platform
import time

import numba
import necklace_enum as te

ROWS = [
    ("fcc binary, all compositions, n=16", "fcc", 16, 2, None),
    ("fcc binary, all compositions, n=20", "fcc", 20, 2, None),
    ("fcc binary, all compositions, n=24", "fcc", 24, 2, None),
    ("fcc ternary, all compositions, n=12", "fcc", 12, 3, None),
    ("fcc ternary, all compositions, n=14", "fcc", 14, 3, None),
    ("fcc quaternary, all compositions, n=10", "fcc", 10, 4, None),
    ("fcc 15:5, n=20", "fcc", 20, 2, (15, 5)),
    ("fcc 24:8, n=32", "fcc", 32, 2, (24, 8)),
    ("fcc 10:5:5, n=20", "fcc", 20, 3, (10, 5, 5)),
    ("hcp binary, all compositions, 24 sites", "hcp", 12, 2, None),
    ("hcp ternary, all compositions, 16 sites", "hcp", 8, 3, None),
]


def timed(f, cooldown):
    time.sleep(cooldown)
    t = time.perf_counter()
    r = f()
    return r, time.perf_counter() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cooldown", type=float, default=0,
                    help="seconds to pause before each timed run, so a fanless machine starts every run cool")
    ap.add_argument("--quick", action="store_true", help="only the first and fourth problems")
    a = ap.parse_args()
    rows = [ROWS[0], ROWS[3]] if a.quick else ROWS
    print(f"# {platform.platform()}, Python {platform.python_version()}, necklace-enum {te.__version__}, "
          f"Numba {numba.__version__}, cooldown {a.cooldown:g} s")
    print(f"{'problem':42s} {'structures':>12s} {'count':>8s} {'arrays':>8s} {'objects':>8s}")
    for name, parent, n, k, comp in rows:
        P = getattr(te, parent)(k)
        kw = dict(composition=comp, label_exchange=False)
        count = lambda: te.count_structures(P, n, **kw)[n]
        arrays = lambda: sum(lab.shape[0] for _, lab in te.enumerate_labelings(P, n, copy=True, **kw))
        objects = lambda: sum(1 for _ in te.enumerate_structures(P, n, degeneracy=False, **kw))
        count()                                                  # warm-up / compilation
        c, tc = None, None
        for _ in range(3):
            c, t = timed(count, a.cooldown)
            tc = t if tc is None else min(tc, t)
        na, ta = timed(arrays, a.cooldown)
        no, to = timed(objects, a.cooldown)
        flag = "" if na == no == c else f"  MISMATCH: arrays {na:,d}, objects {no:,d}"
        print(f"{name:42s} {c:>12,d} {tc:8.2f} {ta:8.2f} {to:8.2f}{flag}", flush=True)


if __name__ == "__main__":
    main()
