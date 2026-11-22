"""Time enumlib (enum.x) and necklace-enum on the problems of the paper's enumlib table, on this machine.

    python3 scripts/time_enumlib.py ~/Projects/enumlib/src/enum.x            # all rows
    python3 scripts/time_enumlib.py ~/Projects/enumlib/src/enum.x --timeout 3600

Both programs do the whole job: superlattices, symmetry, labelings. enumlib's output file is a
link to /dev/null, so it still formats every structure but writes nothing to disk; necklace-enum
counts without listing. enumlib times are single wall-clock runs, with peak memory from
/usr/bin/time; necklace-enum times are the best of three runs after a warm-up run (Numba
compilation excluded). No row uses label exchange ("full list" in enumlib), the same problems and
conventions as the dsenum and Enumlib.jl comparisons, so all four codes go in one table. The ratio
compares enumlib with necklace-enum's listing time (every structure handed to Python).
"""
import argparse
import os
import platform
import re
import subprocess
import sys
import tempfile
import time

import necklace_enum as te

C = (8 / 3) ** 0.5
SQ3 = 3 ** 0.5

# name, parent, cells n, species k, composition (atoms per cell) or None
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


def enum_input(parent, n, k, comp):
    lab = "/".join(map(str, range(k)))
    mode = "full"                          # no label exchange, as in the dsenum and Enumlib.jl runs
    if parent == "fcc":
        s = (f"fcc\nbulk\n0.5 0.5 0.0\n0.5 0.0 0.5\n0.0 0.5 0.5\n {k} -nary case\n 1\n 0.0 0.0 0.0 {lab}\n"
             f" {n} {n}\n0.10000000E-06\n{mode} list of labelings\n")
        N = n
    else:
        s = (f"hcp\nbulk\n1 0 0\n0.5 {SQ3 / 2:.12f} 0\n0 0 {C:.12f}\n {k} -nary case\n 2\n"
             f" 0 0 0 {lab}\n 0.5 {SQ3 / 6:.12f} {C / 2:.12f} {lab}\n"
             f" {n} {n}\n0.0001\n{mode} list of labelings\n")
        N = 2 * n
    if comp:
        s += "# Concentration restrictions\n" + "\n".join(f"{c} {c} {N}" for c in comp) + "\n"
    return s


def run_enumlib(enumx, parent, n, k, comp, timeout):
    work = tempfile.mkdtemp(prefix="enumlib_")
    open(os.path.join(work, "struct_enum.in"), "w").write(enum_input(parent, n, k, comp))
    os.symlink("/dev/null", os.path.join(work, "struct_enum.out"))
    timer = ["/usr/bin/time", "-l"] if platform.system() == "Darwin" else ["/usr/bin/time", "-v"]
    if not os.path.exists("/usr/bin/time"):
        timer = []
    t = time.perf_counter()
    try:
        p = subprocess.run(timer + [enumx], cwd=work, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, None, "timeout"
    wall = time.perf_counter() - t
    m = re.search(r"(\d+)\s+peak memory footprint", p.stderr) or \
        re.search(r"(\d+)\s+maximum resident set size", p.stderr) or \
        re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", p.stderr)
    mem = None
    if m:
        mem = int(m.group(1)) / (1e9 if platform.system() == "Darwin" else 1e6)   # GB
    status = "ok" if p.returncode == 0 else f"exit {p.returncode}"
    return wall, mem, status


def run_necklace(parent, n, k, comp):
    """Count (best of three) and list (one run) without label exchange; warm-up first."""
    P = getattr(te, parent)(k)
    f = lambda: te.count_structures(P, n, composition=comp, label_exchange=False)[n]
    f()
    best, c = None, None
    for _ in range(3):
        t = time.perf_counter(); c = f(); dt = time.perf_counter() - t
        best = dt if best is None else min(best, dt)
    t = time.perf_counter()
    sum(1 for _ in te.enumerate_structures(P, n, composition=comp, label_exchange=False))
    return c, best, time.perf_counter() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("enumx")
    ap.add_argument("--timeout", type=float, default=1800)
    a = ap.parse_args()
    enumx = os.path.abspath(os.path.expanduser(a.enumx))
    if not os.access(enumx, os.X_OK):
        sys.exit(f"not an executable: {enumx}")
    te.count_structures(te.fcc(2), 4)
    print(f"# {platform.platform()}, Python {platform.python_version()}")
    print(f"{'problem':42s} {'structures':>12s} {'enumlib (s)':>12s} {'mem (GB)':>9s} {'count (s)':>10s} {'list (s)':>9s} {'vs list':>8s}")
    for name, parent, n, k, comp in ROWS:
        c, tc, tn = run_necklace(parent, n, k, comp)
        te_, mem, status = run_enumlib(enumx, parent, n, k, comp, a.timeout)
        el = f"{te_:12.2f}" if te_ is not None else f"{status:>12s}"
        if status not in ("ok", "timeout"):
            el = f"{status:>12s}"
        ms = f"{mem:9.2f}" if mem is not None else f"{'-':>9s}"
        ratio = f"{'-':>7s}"
        if te_ is not None and status == "ok":
            r = te_ / tn
            ratio = f"{r:7.1f}" if r < 10 else f"{r:7.0f}"
        print(f"{name:42s} {c:>12,d} {el} {ms} {tc:10.2f} {tn:9.2f} {ratio}", flush=True)


if __name__ == "__main__":
    main()
