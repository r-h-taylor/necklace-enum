"""Time count_structures with and without inversion pruning (Table 1 problems)."""
import time
import necklace_enum as te
import necklace_enum.core as core

ROWS = [("fcc binary, all compositions, n=24", te.fcc(2), 24, None, True),
        ("fcc binary, all compositions, n=26", te.fcc(2), 26, None, True),
        ("fcc ternary, all compositions, n=16", te.fcc(3), 16, None, True),
        ("fcc quaternary, all compositions, n=12", te.fcc(4), 12, None, True),
        ("fcc 10:5:5, n=20", te.fcc(3), 20, (10, 5, 5), False),
        ("fcc 24:8, n=32", te.fcc(2), 32, (24, 8), False),
        ("hcp ternary, all compositions, 16 sites", te.hcp(3), 8, None, True),
        ("hcp 12:12, 24 sites", te.hcp(2), 12, (12, 12), False)]

def best(f, r=3):
    f(); ts = []
    for _ in range(r):
        t = time.perf_counter(); out = f(); ts.append(time.perf_counter() - t)
    return out, min(ts)

on = core.split_reversal
off = lambda uni, members, m, a: (None, uni, members)
print(f"{'problem':42s} {'with (s)':>9s} {'without (s)':>11s} {'factor':>7s}")
for name, P, n, comp, le in ROWS:
    f = lambda: te.count_structures(P, n, composition=comp, label_exchange=le)[n]
    core.split_reversal = on;  c1, t1 = best(f)
    core.split_reversal = off; c2, t2 = best(f)
    core.split_reversal = on
    assert c1 == c2, (name, c1, c2)
    print(f"{name:42s} {t1:9.2f} {t2:11.2f} {t2/t1:7.2f}", flush=True)
