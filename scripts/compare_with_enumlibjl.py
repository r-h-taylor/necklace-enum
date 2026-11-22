"""Compare necklace-enum with Enumlib.jl (Hart's Julia reimplementation of enumlib) on the
problems of the dsenum table.

    python3 scripts/compare_with_enumlibjl.py
    python3 scripts/compare_with_enumlibjl.py --julia ~/.juliaup/bin/julia --timeout 3600
    python3 scripts/compare_with_enumlibjl.py --skip-resource-check

Enumlib.jl has no label-exchange option, so both codes run without it. Enumlib.jl returns every
structure, so its time is compared with necklace-enum's listing time (every structure handed to
Python); counting times are shown too. Each Enumlib.jl run is a separate Julia process, timed
inside Julia after a warm-up call (compilation and package loading excluded). Enumlib.jl's
default resource check may refuse very large runs; these are reported as "refused".
Whether structures lacking an atom type are counted is matched to Enumlib.jl's count.
"""
import argparse
import os
import platform
import subprocess
import time

import necklace_enum as te

HERE = os.path.dirname(os.path.abspath(__file__))
JL = os.path.join(HERE, "time_enumlibjl.jl")

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


def best(f, repeat):
    f()
    t_best, out = None, None
    for _ in range(repeat):
        t = time.perf_counter(); out = f(); dt = time.perf_counter() - t
        t_best = dt if t_best is None else min(t_best, dt)
    return out, t_best


def enumlibjl(julia, parent, n, k, comp, timeout, skip):
    args = [julia, "--startup-file=no", JL, parent, str(n), str(k), "all" if comp is None else ",".join(map(str, comp))]
    if skip:
        args.append("--skip-resource-check")
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, None, "timeout"
    for line in p.stdout.splitlines():
        if line.startswith("RESULT"):
            _, c, t = line.split()
            return int(c), float(t), "ok"
        if line.startswith("REFUSED"):
            return None, None, "refused"
        if line.startswith("ERROR"):
            return None, None, "error: " + line[6:80]
    return None, None, "failed: " + (p.stderr.strip().splitlines() or ["?"])[-1][:80]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--julia", default="julia")
    ap.add_argument("--timeout", type=float, default=1800)
    ap.add_argument("--skip-resource-check", action="store_true")
    a = ap.parse_args()
    v = subprocess.run([a.julia, "--startup-file=no", "-e", "using Enumlib; print(pkgversion(Enumlib), \" Julia \", VERSION)"],
                       capture_output=True, text=True).stdout.strip()
    te.count_structures(te.fcc(2), 4)
    print(f"# {platform.platform()}, Python {platform.python_version()}, Enumlib.jl {v}")
    print(f"{'problem':42s} {'structures':>12s} {'Enumlib.jl':>11s} {'count':>7s} {'list':>7s} {'vs list':>8s}")
    for name, parent, n, k, comp in ROWS:
        cj, tj, status = enumlibjl(a.julia, parent, n, k, comp, a.timeout, a.skip_resource_check)
        P = getattr(te, parent)(k)
        inc = False
        if cj is not None and k > 2 and comp is None:
            c0 = te.count_structures(P, n, label_exchange=False)[n]
            inc = cj != c0                       # Enumlib.jl may include structures lacking a type
        cnt, tc = best(lambda: te.count_structures(P, n, composition=comp, label_exchange=False,
                                                   include_incomplete=inc)[n], 3)
        _, tl = best(lambda: sum(1 for _ in te.enumerate_structures(P, n, composition=comp, label_exchange=False,
                                                                     include_incomplete=inc)), 1)
        flag = "" if (cj is None or cj == cnt) else f"  COUNT MISMATCH: Enumlib.jl {cj:,d}"
        tjs = f"{tj:11.2f}" if tj is not None else f"{status:>11s}"
        ratio = f"{tj / tl:8.1f}" if tj is not None else f"{'-':>8s}"
        note = " (incl. structures lacking a type)" if inc else ""
        print(f"{name:42s} {cnt:>12,d} {tjs} {tc:7.2f} {tl:7.2f} {ratio}{flag}{note}", flush=True)


if __name__ == "__main__":
    main()
