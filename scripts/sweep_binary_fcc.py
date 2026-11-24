"""Size sweep for the scaling figure: binary fcc, all compositions, no label exchange.

    caffeinate -i python3 scripts/sweep_binary_fcc.py --cooldown 60 > scripts/sweep_binary_fcc_M2.txt

Every run is a separate process under /usr/bin/time -l, which gives its peak memory (macOS peak
memory footprint). The listing time is measured inside the process, as in Table 1:
  necklace   enumerate_labelings (structures handed to Python as integer arrays), after a warm-up
             on a small cell so that Numba compilation is excluded
  dsenum     building the decision diagram and listing every labeling (no warm-up, as in Table 1)
  enumlibjl  enumerate_structures, timed inside Julia after a warm-up on a small cell
  enumlib    wall-clock time of enum.x, writing its formatted output to /dev/null
All codes run n = 10, 12, ..., 24; the present method also runs n = 26 and 28.
With --composition 3,1 every run is at that fixed composition instead (atoms per cell scaled to n):

    caffeinate -i python3 scripts/sweep_binary_fcc.py --composition 3,1 --sizes 16,20,24,28,32 \
        --extra 36 --cooldown 60 > scripts/sweep_binary_fcc_31_M2.txt
Each run stops after --timeout seconds; --cooldown pauses before every run.
"""
import argparse
import os
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from memory_compare import timed, enumlib_run, julia_cmd  # noqa: E402

NECKLACE = """
import time, necklace_enum as te
P = te.fcc(2)
sum(l.shape[0] for _, l in te.enumerate_labelings(P, 4, label_exchange=False))   # warm-up
t = time.perf_counter()
c = sum(l.shape[0] for _, l in te.enumerate_labelings(P, {n}, composition={comp!r}, label_exchange=False))
print("RESULT", c, time.perf_counter() - t)
"""

DSENUM = """
import sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, {here!r})
from compare_with_dsenum import ds_list, fcc_base
B = fcc_base()
t = time.perf_counter()
c = ds_list(B, {n}, 2, {comp!r})
print("RESULT", c, time.perf_counter() - t)
"""


def parse(out):
    for line in out.splitlines():
        if line.startswith("RESULT"):
            _, c, t = line.split()
            return int(c), float(t)
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enumx", default=os.path.expanduser("~/Projects/enumlib/src/enum.x"))
    ap.add_argument("--julia", default="julia")
    ap.add_argument("--codes", default="necklace,dsenum,enumlibjl,enumlib")
    ap.add_argument("--sizes", default="10,12,14,16,18,20,22,24")
    ap.add_argument("--extra", default="26,28", help="sizes run by the present method only")
    ap.add_argument("--composition", default=None, help="fixed composition as a ratio, e.g. 3,1")
    ap.add_argument("--timeout", type=float, default=3600)
    ap.add_argument("--cooldown", type=float, default=0)
    a = ap.parse_args()
    codes = a.codes.split(",")
    sizes = [int(x) for x in a.sizes.split(",")]
    extra = [int(x) for x in a.extra.split(",")] if a.extra else []
    if "enumlibjl" in a.codes:
        import subprocess as _sp
        _v = _sp.run([a.julia, "--startup-file=no", "-e", 'using Enumlib; print(pkgversion(Enumlib), " Julia ", VERSION)'],
                     capture_output=True, text=True, stdin=_sp.DEVNULL).stdout.strip()
        print(f"# Enumlib.jl {_v}, JULIA_PROJECT={os.environ.get('JULIA_PROJECT', '(default)')}", flush=True)
    print(f"# {platform.platform()}, Python {platform.python_version()}, cooldown {a.cooldown:g} s, "
          f"timeout {a.timeout:g} s", flush=True)
    print(f"{'code':10s} {'n':>3s} {'structures':>12s} {'time_s':>10s} {'peak_GB':>8s}  status", flush=True)
    jobs = [(n, c) for n in sizes for c in codes] + [(n, "necklace") for n in extra if "necklace" in codes]
    ratio = [int(x) for x in a.composition.split(",")] if a.composition else None
    if ratio:
        print(f"# fixed composition {a.composition} (atoms per cell scaled to n)", flush=True)
    for n, code in jobs:
        comp = None
        if ratio:
            if n % sum(ratio):
                continue
            comp = tuple(r * n // sum(ratio) for r in ratio)
        time.sleep(a.cooldown)
        if code == "necklace":
            mem, wall, status, out = timed([sys.executable, "-c", NECKLACE.format(n=n, comp=comp)], timeout=a.timeout)
            c, t = parse(out)
        elif code == "dsenum":
            mem, wall, status, out = timed([sys.executable, "-c", DSENUM.format(n=n, here=HERE, comp=comp)], timeout=a.timeout)
            c, t = parse(out)
        elif code == "enumlibjl":
            mem, wall, status, out = timed(julia_cmd(a.julia, "fcc", n, 2, comp), timeout=a.timeout)
            c, t = parse(out)
        else:
            mem, wall, status, out = enumlib_run(a.enumx, "fcc", n, 2, comp, a.timeout)
            c, t = None, wall
        cs = f"{c:,d}" if c is not None else "-"
        ts = f"{t:.3f}" if t is not None else "-"
        ms = f"{mem:.3f}" if mem is not None else "-"
        print(f"{code:10s} {n:3d} {cs:>12s} {ts:>10s} {ms:>8s}  {status}", flush=True)


if __name__ == "__main__":
    main()
