"""Command line: count or list derivative structures.

Examples
--------
    necklace-enum count fcc Cu,Au 1-12
    necklace-enum count fcc 3 2-10 --label-exchange
    necklace-enum count fcc 2 1-60 --label-exchange --burnside   # counts only, no listing
    necklace-enum list hcp Mg,Cd 4 --composition 5,3 --format extxyz -o out.xyz
    necklace-enum list fcc Cu,Au 1-8 --format enum -o struct_enum.out
    necklace-enum list fcc Cu,Au 8 --composition 6,2 --degeneracy   # with each structure's weight g
    necklace-enum enumlib [struct_enum.in]      # drop-in for enumlib's enum.x
"""
from __future__ import annotations

import argparse
import sys
import time

from . import bcc, count_structures, count_structures_burnside, enumerate_structures, fcc, hcp, sc

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
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "enumlib":
        from .enumlib_format import run_enumlib_input
        infile = argv[1] if len(argv) > 1 else "struct_enum.in"
        t = time.perf_counter()
        n = run_enumlib_input(infile, "struct_enum.out")
        print(f"wrote struct_enum.out: {n:,d} structures ({time.perf_counter() - t:.1f} s)")
        return 0
    ap = argparse.ArgumentParser(prog="necklace-enum", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["count", "list"])
    ap.add_argument("parent", choices=sorted(PARENTS))
    ap.add_argument("species", help="comma-separated names (Cu,Au) or a number of species")
    ap.add_argument("sizes", help="supercell sizes, e.g. 1-12 or 4,6,8")
    ap.add_argument("--composition", help="atom counts per species, e.g. 6,2")
    ap.add_argument("--label-exchange", action="store_true",
                    help="identify structures related by permuting species (Hart-Forcade tables)")
    ap.add_argument("--format", default="summary",
                    help="'summary', 'enum' (enumlib's struct_enum.out) or an ASE format such as extxyz")
    ap.add_argument("--include-incomplete", action="store_true",
                    help="also list structures in which some species is absent")
    ap.add_argument("--burnside", action="store_true",
                    help="count by Burnside's lemma instead of generating (count only; reaches far larger sizes)")
    ap.add_argument("--degeneracy", action="store_true",
                    help="summary: also print each structure's degeneracy on its superlattice, the number of "
                         "equivalent superlattices and their product (the weight g in a partition function)")
    ap.add_argument("-o", "--output", help="output file for --format other than summary")
    args = ap.parse_args(argv)

    P = PARENTS[args.parent](_species(args.species))
    comp = [int(x) for x in args.composition.split(",")] if args.composition else None
    sizes = _sizes(args.sizes)
    if args.action == "count":
        total = 0
        for n in sizes:
            t = time.perf_counter()
            counter = count_structures_burnside if args.burnside else count_structures
            c = counter(P, n, composition=comp, label_exchange=args.label_exchange,
                        include_incomplete=args.include_incomplete)[n]
            total += c
            print(f"n={n:3d}  {c:>14,d}  ({time.perf_counter() - t:.2f} s)")
        print(f"total {total:,d}")
        return 0
    if args.format == "enum":
        from .enumlib_format import write_struct_enum
        n = write_struct_enum(args.output or "struct_enum.out", P, sizes, composition=comp,
                              include_incomplete=args.include_incomplete)
        print(f"wrote {args.output or 'struct_enum.out'}: {n:,d} structures")
        return 0
    structs = enumerate_structures(P, sizes, composition=comp, label_exchange=args.label_exchange,
                                   include_incomplete=args.include_incomplete)
    if args.format == "summary":
        if args.degeneracy:
            print("# size  HNF  labeling  degeneracy  superlattice_degeneracy  total_degeneracy")
        for s in structs:
            row = [s.size, s.hnf.ravel().tolist(), "".join(str(P.species[x])[0] for x in s.labels)]
            if args.degeneracy:
                row += [s.degeneracy, s.superlattice_degeneracy, s.total_degeneracy]
            print(*row)
        return 0
    try:
        from ase.io import write
    except ImportError:
        print(f"necklace-enum: --format {args.format} writes structures with ASE, which is not installed.\n"
              "Install it with:  pip install ase", file=sys.stderr)
        return 1
    write(args.output or sys.stdout, [s.to_ase() for s in structs], format=args.format)
    return 0


if __name__ == "__main__":
    sys.exit(main())
