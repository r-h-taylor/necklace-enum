"""Write OEIS b-files of derivative-structure counts, using the Burnside counter.

    python scripts/oeis_counts.py                        # binary fcc/bcc, up to label exchange, n = 1..100
    python scripts/oeis_counts.py --nmax 60 --no-label-exchange
    python scripts/oeis_counts.py --parent sc --species 3 --nmax 40

Each line of the b-file is "n a(n)". Terms up to n = 32 can be checked against the
generator (count_structures); see tests/test_burnside.py. The timing of each term
goes to stderr.
"""
import argparse
import sys
import time

import necklace_enum as te


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parent", default="fcc", choices=["fcc", "bcc", "sc", "hcp"])
    ap.add_argument("--species", type=int, default=2)
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=100)
    ap.add_argument("--no-label-exchange", action="store_true")
    ap.add_argument("-o", "--output", help="b-file name (default: stdout)")
    a = ap.parse_args()
    P = getattr(te, a.parent)(a.species)
    out = open(a.output, "w") if a.output else sys.stdout
    for n in range(a.nmin, a.nmax + 1):
        t = time.perf_counter()
        c = te.count_structures_burnside(P, n, label_exchange=not a.no_label_exchange)[n]
        print(n, c, file=out, flush=True)
        print(f"n={n}: {time.perf_counter() - t:.1f} s", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
