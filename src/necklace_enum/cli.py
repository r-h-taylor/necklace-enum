"""Command line: count or list derivative structures.

Examples
--------
    necklace-enum count fcc Cu,Au 1-12
    necklace-enum count fcc 3 2-10 --label-exchange
    necklace-enum list hcp Mg,Cd 4 --composition 5,3 --format extxyz -o out.xyz
"""
from __future__ import annotations

import argparse
import sys
import time

from . import bcc, count_structures, enumerate_structures, fcc, hcp, sc

PARENTS = {"fcc": fcc, "bcc": bcc, "sc": sc, "hcp": hcp}


def _sizes(text):
    out = []
    for part in text.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def _species(text):
    return int(text) if text.isdigit() else text.split(",")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="necklace-enum", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["count", "list"])
    ap.add_argument("parent", choices=sorted(PARENTS))
    ap.add_argument("species", help="comma-separated names (Cu,Au) or a number of species")
    ap.add_argument("sizes", help="supercell sizes, e.g. 1-12 or 4,6,8")
    ap.add_argument("--composition", help="atom counts per species, e.g. 6,2")
    ap.add_argument("--label-exchange", action="store_true",
                    help="identify structures related by permuting species (Hart-Forcade tables)")
    ap.add_argument("--format", default="summary", help="'summary' or an ASE format such as extxyz")
    ap.add_argument("-o", "--output", help="output file for --format other than summary")
    args = ap.parse_args(argv)

    P = PARENTS[args.parent](_species(args.species))
    comp = [int(x) for x in args.composition.split(",")] if args.composition else None
    sizes = _sizes(args.sizes)
    if args.action == "count":
        total = 0
        for n in sizes:
            t = time.perf_counter()
            c = count_structures(P, n, composition=comp, label_exchange=args.label_exchange)[n]
            total += c
            print(f"n={n:3d}  {c:>14,d}  ({time.perf_counter() - t:.2f} s)")
        print(f"total {total:,d}")
        return 0
    structs = enumerate_structures(P, sizes, composition=comp, label_exchange=args.label_exchange)
    if args.format == "summary":
        for s in structs:
            print(s.size, s.hnf.ravel().tolist(), "".join(str(P.species[x])[0] for x in s.labels))
        return 0
    from ase.io import write
    write(args.output or sys.stdout, [s.to_ase() for s in structs], format=args.format)
    return 0


if __name__ == "__main__":
    sys.exit(main())
