"""Run enumlib's enum.x and necklace-enum on the same struct_enum.in and compare.

    python scripts/compare_with_enumlib.py /path/to/enum.x struct_enum.in

Writes enumlib/struct_enum.out and necklace-enum/struct_enum.out next to the input
and checks, for every cell size, that the number of structures and the multiset of
superlattice data (HNF degeneracy, number of fixing operations, composition up to
renaming of species) agree. The two programs may choose different but equivalent
representatives, so rows are not compared line by line. enumlib's labdegn column is
not compared: it is 0 in its fixed-concentration mode and not an orbit size for some
labelings in its "full" mode.
"""
import collections
import os
import shutil
import subprocess
import sys
import time

import necklace_enum as te


def signature(rows):
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        comp = tuple(sorted(collections.Counter(r["labeling"]).values()))
        by[r["n"]][(r["hnf_degen"], r["pg"], comp)] += 1
    return by


def main(enumx, infile):
    base = os.path.dirname(os.path.abspath(infile))
    dirs = {name: os.path.join(base, name) for name in ("enumlib", "necklace-enum")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
        shutil.copy(infile, os.path.join(d, "struct_enum.in"))
    t = time.perf_counter()
    subprocess.run([os.path.abspath(enumx)], cwd=dirs["enumlib"], check=True, stdout=subprocess.DEVNULL)
    t_el = time.perf_counter() - t
    t = time.perf_counter()
    te.run_enumlib_input(os.path.join(dirs["necklace-enum"], "struct_enum.in"),
                         os.path.join(dirs["necklace-enum"], "struct_enum.out"))
    t_ne = time.perf_counter() - t
    _, a = te.read_struct_enum(os.path.join(dirs["necklace-enum"], "struct_enum.out"))
    _, b = te.read_struct_enum(os.path.join(dirs["enumlib"], "struct_enum.out"))
    sa, sb = signature(a), signature(b)
    print(f"necklace-enum: {len(a):,d} structures in {t_ne:.2f} s")
    print(f"enumlib:       {len(b):,d} structures in {t_el:.2f} s")
    ok = True
    for n in sorted(set(sa) | set(sb)):
        same = sa[n] == sb[n]
        ok &= same
        print(f"  n={n:3d}: {sum(sa[n].values()):8d} vs {sum(sb[n].values()):8d}  {'same' if same else 'DIFFERENT'}")
    print("outputs agree" if ok else "OUTPUTS DIFFER")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
