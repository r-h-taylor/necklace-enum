"""Peak memory (on macOS the peak memory footprint, which includes compressed and swapped memory) of necklace-enum, dsenum, Enumlib.jl and enumlib on the
problems of the paper's Table 1, each code in its own process, measured with /usr/bin/time -l.

    python3 scripts/memory_compare.py --enumx ~/Projects/enumlib/src/enum.x
    python3 scripts/memory_compare.py --enumx ~/Projects/enumlib/src/enum.x --big     # adds n=32 and 48-site rows
    python3 scripts/memory_compare.py --codes necklace,dsenum                         # a subset of codes
    python3 scripts/memory_compare.py --codes dsenum,enumlibjl --only-big --timeout 3600

A run that ends in an allocation error, or is killed by the operating system, is reported as OOM,
with its peak memory when /usr/bin/time could still record it; a run stopped at --timeout as timeout.

For each code a baseline run on a tiny cell is reported first (interpreter, libraries and
compiled code); "growth" in the analysis is peak minus that baseline. necklace-enum is measured
both counting and listing (every structure consumed as it is generated, nothing kept).
No label exchange, as in Table 1. Each run stops after --timeout seconds.
"""
import argparse
import os
import platform
import re
import signal
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
JL = os.path.join(HERE, "time_enumlibjl.jl")
C = (8 / 3) ** 0.5
SQ3 = 3 ** 0.5

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
BIG = [
    ("fcc binary, all compositions, n=28", "fcc", 28, 2, None),
    ("fcc binary, all compositions, n=32", "fcc", 32, 2, None),
    ("fcc 44:4, n=48", "fcc", 48, 2, (44, 4)),
    ("fcc 42:6, n=48", "fcc", 48, 2, (42, 6)),
]
BASE = ("baseline (tiny cell)", "fcc", 2, 2, None)



def run_group(cmd, cwd=None, timeout=None):
    """subprocess.run with a timeout that also kills grandchildren (e.g. enum.x under /usr/bin/time)."""
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                         start_new_session=True)
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        p.communicate()
        raise
    return subprocess.CompletedProcess(cmd, p.returncode, out, err)

def timed(cmd, cwd=None, timeout=1800):
    """Run cmd under /usr/bin/time; return (peak GB, wall s, status, stdout)."""
    darwin = platform.system() == "Darwin"
    full = (["/usr/bin/time", "-l"] if darwin else ["/usr/bin/time", "-v"]) + cmd
    try:
        p = run_group(full, cwd=cwd, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, None, "timeout", ""
    # macOS: "peak memory footprint" counts memory the system has compressed or swapped out,
    # which the resident set size does not; it is the better measure when a run exceeds RAM.
    m = re.search(r"(\d+)\s+peak memory footprint", p.stderr) or \
        re.search(r"(\d+)\s+maximum resident set size", p.stderr) or \
        re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", p.stderr)
    w = re.search(r"([\d.]+) real", p.stderr) or re.search(r"Elapsed \(wall clock\) time.*: ([\d:.]+)", p.stderr)
    mem = int(m.group(1)) / (1e9 if darwin else 1e6) if m else None
    wall = float(w.group(1)) if (w and ":" not in w.group(1)) else None
    out = p.stdout + p.stderr
    failed = p.returncode != 0 or re.search(r"^(ERROR|REFUSED)", p.stdout, re.M)
    # out of memory: an allocation error from the program itself, or a kill by the operating
    # system (SIGKILL that is not ours: our own time limit is reported as "timeout" above)
    oom = re.search(r"OutOfMemory|MemoryError|out of memory|Cannot allocate memory|bad_alloc", out, re.I) \
        or p.returncode in (-9, 137) or re.search(r"signal 9", p.stderr)
    if not failed:
        status = "ok"
    elif oom:
        status = "out of memory"
    else:
        last = (p.stdout.strip() or p.stderr.strip()).splitlines()
        status = f"exit {p.returncode}: {last[-1][:60] if last else ''}"
    return mem, wall, status, p.stdout


def py(code):
    return [sys.executable, "-c", code]


def necklace_cmd(parent, n, k, comp, mode):
    act = (f"print(te.count_structures(P, {n}, composition={comp!r}, label_exchange=False)[{n}])" if mode == "count" else
           f"print(sum(1 for _ in te.enumerate_structures(P, {n}, composition={comp!r}, label_exchange=False)))")
    return py(f"import necklace_enum as te; P = te.{parent}({k}); {act}")


def dsenum_cmd(parent, n, k, comp):
    base = ("Structure(Lattice([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]), ['Cu'], [[0,0,0]])" if parent == "fcc" else
            "Structure(Lattice.hexagonal(1.0, (8/3)**0.5), ['Mg','Mg'], [[1/3,2/3,0.25],[2/3,1/3,0.75]])")
    ratio = None
    if comp:
        from math import gcd
        from functools import reduce
        g = reduce(gcd, comp); ratio = [c // g for c in comp]
    return py("import warnings; warnings.filterwarnings('ignore'); "
              "from pymatgen.core import Structure, Lattice; from dsenum import ZddStructureEnumerator; "
              f"e = ZddStructureEnumerator({base}, {n}, {k}, composition_constraints={ratio!r}, "
              "remove_superperiodic=True, remove_incomplete=True, verbose=False); print(e.count())")


def julia_cmd(julia, parent, n, k, comp):
    return [julia, "--startup-file=no", JL, parent, str(n), str(k), "all" if comp is None else ",".join(map(str, comp)),
            "--skip-resource-check"]


def enumlib_run(enumx, parent, n, k, comp, timeout):
    lab = "/".join(map(str, range(k)))
    if parent == "fcc":
        s = (f"fcc\nbulk\n0.5 0.5 0.0\n0.5 0.0 0.5\n0.0 0.5 0.5\n {k} -nary case\n 1\n 0.0 0.0 0.0 {lab}\n"
             f" {n} {n}\n0.10000000E-06\nfull list of labelings\n"); N = n
    else:
        s = (f"hcp\nbulk\n1 0 0\n0.5 {SQ3 / 2:.12f} 0\n0 0 {C:.12f}\n {k} -nary case\n 2\n"
             f" 0 0 0 {lab}\n 0.5 {SQ3 / 6:.12f} {C / 2:.12f} {lab}\n"
             f" {n} {n}\n0.0001\nfull list of labelings\n"); N = 2 * n
    if comp:
        s += "# Concentration restrictions\n" + "\n".join(f"{c} {c} {N}" for c in comp) + "\n"
    w = tempfile.mkdtemp(prefix="enumlib_mem_")
    open(os.path.join(w, "struct_enum.in"), "w").write(s)
    os.symlink("/dev/null", os.path.join(w, "struct_enum.out"))
    return timed([enumx], cwd=w, timeout=timeout)


def fmt(r):
    mem, wall, status, _ = r
    if status == "out of memory":
        return f"{'OOM' if mem is None else f'OOM {mem:.1f}':>14s}"   # peak before the failure, if known
    if status != "ok":
        return f"{status[:13]:>14s}"
    return f"{mem:7.3f} ({wall:5.0f}s)" if wall is not None else f"{mem:7.3f}       "


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enumx", default=os.path.expanduser("~/Projects/enumlib/src/enum.x"))
    ap.add_argument("--julia", default="julia")
    ap.add_argument("--codes", default="necklace,dsenum,enumlibjl,enumlib")
    ap.add_argument("--big", action="store_true")
    ap.add_argument("--only-big", action="store_true", help="only the baseline and the larger problems")
    ap.add_argument("--timeout", type=float, default=1800)
    a = ap.parse_args()
    codes = a.codes.split(",")
    rows = [BASE] + BIG if a.only_big else [BASE] + ROWS + (BIG if a.big else [])
    cols = []
    if "necklace" in codes: cols += ["necklace count", "necklace list"]
    if "dsenum" in codes: cols += ["dsenum count"]
    if "enumlibjl" in codes: cols += ["Enumlib.jl"]
    if "enumlib" in codes: cols += ["enumlib"]
    if "enumlibjl" in a.codes:
        import subprocess as _sp
        _v = _sp.run([a.julia, "--startup-file=no", "-e", 'using Enumlib; print(pkgversion(Enumlib), " Julia ", VERSION)'],
                     capture_output=True, text=True, stdin=_sp.DEVNULL).stdout.strip()
        print(f"# Enumlib.jl {_v}, JULIA_PROJECT={os.environ.get('JULIA_PROJECT', '(default)')}", flush=True)
    print(f"# {platform.platform()}, Python {platform.python_version()}; peak resident memory in GB (wall time)")
    print(f"{'problem':42s} " + " ".join(f"{c:>14s}" for c in cols))
    for name, parent, n, k, comp in rows:
        out = []
        if "necklace" in codes:
            out.append(timed(necklace_cmd(parent, n, k, comp, "count"), timeout=a.timeout))
            out.append(timed(necklace_cmd(parent, n, k, comp, "list"), timeout=a.timeout))
        if "dsenum" in codes:
            out.append(timed(dsenum_cmd(parent, n, k, comp), timeout=a.timeout))
        if "enumlibjl" in codes:
            out.append(timed(julia_cmd(a.julia, parent, n, k, comp), timeout=a.timeout))
        if "enumlib" in codes:
            out.append(enumlib_run(a.enumx, parent, n, k, comp, a.timeout))
        print(f"{name:42s} " + " ".join(fmt(r) for r in out), flush=True)


if __name__ == "__main__":
    main()
