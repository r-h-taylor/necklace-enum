"""Build the paper's Table 1 (LaTeX rows) and the speed ratios quoted in the text from the timing logs.

    python3 scripts/build_table1.py            # reads scripts/*_timings_M2.txt

Inputs (all measured on the same machine with --cooldown):
  enumlib_timings_M2.txt         scripts/time_enumlib.py
  enumlibjl_v050_timings_M2.txt  scripts/compare_with_enumlibjl.py
  dsenum_timings_M2.txt          scripts/compare_with_dsenum.py
  necklace_timings_M2.txt        scripts/time_necklace.py   (count, arrays, objects)
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))

ROWS = [
    ("fcc binary, all compositions, n=16", r"fcc binary, all comp., $n=16$", False),
    ("fcc binary, all compositions, n=20", r"fcc binary, all comp., $n=20$", False),
    ("fcc binary, all compositions, n=24", r"fcc binary, all comp., $n=24$", False),
    ("fcc ternary, all compositions, n=12", r"fcc ternary, all comp., $n=12$", True),
    ("fcc ternary, all compositions, n=14", r"fcc ternary, all comp., $n=14$", True),
    ("fcc quaternary, all compositions, n=10", r"fcc quaternary, all comp., $n=10$", True),
    ("fcc 15:5, n=20", r"fcc 15:5, $n=20$", False),
    ("fcc 24:8, n=32", r"fcc 24:8, $n=32$", False),
    ("fcc 10:5:5, n=20", r"fcc 10:5:5, $n=20$", False),
    ("hcp binary, all compositions, 24 sites", r"hcp binary, all comp., 24 sites", False),
    ("hcp ternary, all compositions, 16 sites", r"hcp ternary, all comp., 16 sites", True),
]


def read(name):
    out = {}
    for line in open(os.path.join(HERE, name)):
        for key, _, _ in ROWS:
            if line.startswith(key + " "):
                out[key] = line[len(key):].split()
    return out


def num(tok):
    try:
        return float(tok.replace(",", ""))
    except ValueError:
        return None                                   # "timeout", "-"


def fmt(x):
    if x is None:
        return r"$>$3600"
    if x < 10:
        return f"{x:.2f}"
    if x < 1000:
        return f"{x:.1f}"
    return f"{x:.0f}"


def main():
    el, jl, ds, ne = (read(f) for f in ("enumlib_timings_M2.txt", "enumlibjl_v050_timings_M2.txt",
                                        "dsenum_timings_M2.txt", "necklace_timings_M2.txt"))
    ratios = {k: [] for k in ("el_obj", "el_arr", "jl_arr", "jl_arr_same", "ds_cnt", "ds_arr", "ds_obj")}
    lines = []
    for key, tex, dagger in ROWS:
        n_struct = ne[key][0]
        t_el = num(el[key][1])
        t_jl = num(jl[key][1])
        t_dsc, t_dsl = num(ds[key][1]), num(ds[key][3])
        t_c, t_a, t_o = (num(x) for x in ne[key][1:4])
        d = r"$^\dagger$" if dagger else ""
        best = min(t for t in (t_el, t_jl, t_dsl) if t is not None)
        sp = best / t_a
        ratios.setdefault("best_arr", []).append((sp, key))
        sps = f"{sp:.1f}" if sp < 10 else f"{sp:.0f}"
        lines.append(f"{tex} & {n_struct.replace(',', chr(92) + ',')} & {fmt(t_el)}{d} & {fmt(t_jl)}{d} & "
                     f"{fmt(t_dsc)} & {fmt(t_dsl)} & {fmt(t_c)} & {fmt(t_a)} & {fmt(t_o)} & {sps}\\\\")
        ratios["el_obj"].append((t_el / t_o, key))
        ratios["el_arr"].append((t_el / t_a, key))
        if t_jl is not None:
            ratios["jl_arr"].append((t_jl / t_a, key))
            if not dagger:
                ratios["jl_arr_same"].append((t_jl / t_a, key))
        ratios["ds_cnt"].append((t_dsc / t_c, key))
        ratios["ds_arr"].append((t_dsl / t_a, key))
        ratios["ds_obj"].append((t_dsl / t_o, key))
    print("\n".join(lines))
    print()
    for k, v in ratios.items():
        if not v:
            continue
        lo, hi = min(v), max(v)
        print(f"{k:12s} {lo[0]:8.2f} ({lo[1]})  ..  {hi[0]:8.1f} ({hi[1]})")
    fcc = [r for r in ratios["ds_obj"] if r[1].startswith("fcc")]
    hcp = [r for r in ratios["ds_obj"] if r[1].startswith("hcp")]
    print(f"ds_obj fcc    {min(fcc)[0]:.1f} .. {max(fcc)[0]:.1f};  hcp {min(hcp)[0]:.2f} .. {max(hcp)[0]:.2f}")
    print(f"count range   {min(num(ne[k][1]) for k, _, _ in ROWS):.2f} .. {max(num(ne[k][1]) for k, _, _ in ROWS):.2f} s")


if __name__ == "__main__":
    main()
