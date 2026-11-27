#!/usr/bin/env python3
"""Check enumlib's labdegn column against an orbit size computed from scratch.

    python3 scripts/check_enumlib_labdegn.py ~/Projects/enumlib/src/enum.x
    python3 scripts/check_enumlib_labdegn.py ~/Projects/enumlib/src/enum.x --keep out_dir

For each case below, enum.x runs on a generated struct_enum.in. Every row of its
struct_enum.out is then read back with enumlib's own conventions (site at lattice
point z has address g = L z mod D, label index g1*d2*d3 + g2*d3 + g3) and the orbit
of that labeling on that superlattice is generated explicitly: every translation of
the supercell, times every parent point operation that maps the superlattice onto
itself, and, as a second number, also every permutation of the species. The orbit
sizes are compared with enumlib's labdegn, and Totdegn with Hdegn * labdegn.

Needs only numpy and enum.x. Shares no code with necklace-enum or enumlib.
Bravais parents only (one site per primitive cell).
"""
import argparse
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

import numpy as np

PARENTS = {
    "fcc": [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],
    "bcc": [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]],
    "sc": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
}

# name, parent, k, sizes (lo, hi), mode, concentration (atom counts at size hi) or None
CASES = [
    ("fcc binary, full", "fcc", 2, (1, 8), "full", None),
    ("fcc binary, part", "fcc", 2, (1, 8), "part", None),
    ("fcc ternary, full", "fcc", 3, (1, 6), "full", None),
    ("fcc ternary, part", "fcc", 3, (1, 6), "part", None),
    ("bcc binary, full", "bcc", 2, (1, 8), "full", None),
    ("sc binary, part", "sc", 2, (1, 8), "part", None),
    ("fcc 6:2 at n=8", "fcc", 2, (8, 8), "full", (6, 2)),
    ("fcc 4:4 at n=8", "fcc", 2, (8, 8), "full", (4, 4)),
    ("fcc 9:3 at n=12", "fcc", 2, (12, 12), "full", (9, 3)),
    ("fcc 2:2:2 at n=6", "fcc", 3, (6, 6), "full", (2, 2, 2)),
]


def point_group(vectors):
    """Integer matrices M (lattice coordinates, columns = images of the basis vectors)
    that preserve the metric A^T A, A with the parent vectors as columns."""
    A = np.array(vectors, float).T
    G = A.T @ A
    out = []
    for e in itertools.product((-1, 0, 1), repeat=9):
        M = np.array(e).reshape(3, 3)
        if round(abs(np.linalg.det(M))) == 1 and np.allclose(M.T @ G @ M, G):
            out.append(M)
    return out


def struct_enum_in(vectors, k, sizes, mode, conc):
    lab = "/".join(str(i) for i in range(k))
    lines = ["labdegn check", "bulk"]
    lines += [" ".join(f"{x:.12f}" for x in v) for v in vectors]
    lines += [f" {k} -nary case", " 1", f" 0.0 0.0 0.0 {lab}", f" {sizes[0]} {sizes[1]}",
              "0.10000000E-06", f"{mode} list of labelings"]
    if conc:
        lines.append("# Concentration restrictions")
        lines += [f"{c} {c} {sizes[1]}" for c in conc]
    return "\n".join(lines) + "\n"


def read_rows(path):
    lines = open(path).read().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("start"))
    rows = []
    for l in lines[start + 1:]:
        f = l.split()
        if len(f) < 27:
            continue
        v = [int(x) for x in f[:26]]
        rows.append({"hnf_degen": v[2], "lab_degen": v[3], "tot_degen": v[4], "n": v[6], "pg": v[7],
                     "D": np.array(v[8:11]),
                     "H": np.array([[v[11], 0, 0], [v[12], v[13], 0], [v[14], v[15], v[16]]]),
                     "L": np.array(v[17:26]).reshape(3, 3), "labeling": f[26], "line": l.strip()})
    return rows


class Cell:
    """The n sites of one superlattice, addressed as enumlib addresses them."""

    def __init__(self, H, L, D, rots):
        self.L, self.D = L, D
        self.n = int(D.prod())
        Hi = np.linalg.inv(H)
        stab = [M for M in rots if np.allclose(Hi @ M @ H, np.round(Hi @ M @ H))]
        rep = {}
        r = int(np.abs(H).max()) + 1
        for z in itertools.product(range(-r, r + 1), repeat=3):
            z = np.array(z)
            rep.setdefault(self.index(z), z)
            if len(rep) == self.n:
                break
        assert len(rep) == self.n, "could not find a lattice point for every address"
        self.rep = [rep[i] for i in range(self.n)]
        # every element of (translations x stabilizer) as a site permutation:
        # the image labeling is lab'(z) = lab(M^-1 (z - t)), so lab'[i] = lab[perm[i]]
        self.perms = set()
        for M in stab:
            Mi = np.round(np.linalg.inv(M)).astype(int)
            for t in self.rep:
                self.perms.add(tuple(self.index(Mi @ (z - t)) for z in self.rep))

    def index(self, z):
        g = (self.L @ z) % self.D
        return int(g[0] * self.D[1] * self.D[2] + g[1] * self.D[2] + g[2])

    def orbit(self, lab, label_perms):
        """Number of distinct labelings equivalent to lab."""
        seen = set()
        for p in self.perms:
            img = [lab[i] for i in p]
            for e in label_perms:
                seen.add(tuple(e[x] for x in img))
        return len(seen)


def run_case(enumx, workdir, name, parent, k, sizes, mode, conc, show):
    os.makedirs(workdir, exist_ok=True)
    vectors = PARENTS[parent]
    with open(os.path.join(workdir, "struct_enum.in"), "w") as fh:
        fh.write(struct_enum_in(vectors, k, sizes, mode, conc))
    r = subprocess.run([enumx], cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    out = os.path.join(workdir, "struct_enum.out")
    if r.returncode != 0 or not os.path.exists(out):
        print(f"{name:22s} enum.x failed (exit {r.returncode}): {r.stderr.strip()[:200]}")
        return None
    rots = point_group(vectors)
    no_ex = [tuple(range(k))]
    all_ex = list(itertools.permutations(range(k)))
    cells = {}
    stats = Counter()
    bad = []
    for row in read_rows(out):
        key = (row["H"].tobytes(), row["L"].tobytes())
        if key not in cells:
            cells[key] = Cell(row["H"], row["L"], row["D"], rots)
        cell = cells[key]
        lab = [int(c) for c in row["labeling"]]
        o1 = cell.orbit(lab, no_ex)
        o2 = cell.orbit(lab, all_ex)
        got = row["lab_degen"]
        stats["rows"] += 1
        stats["zero"] += got == 0
        stats["= orbit"] += got == o1
        stats["= orbit with exchange"] += got == o2
        stats["Totdegn = Hdegn*labdegn"] += row["tot_degen"] == row["hnf_degen"] * got
        if got not in (o1, o2):
            stats["neither"] += 1
            bad.append((row, o1, o2))
    s = stats
    print(f"{name:22s} {s['rows']:6d} {s['zero']:6d} {s['= orbit']:8d} {s['= orbit with exchange']:10d} "
          f"{s['neither']:8d} {s['Totdegn = Hdegn*labdegn']:8d}")
    for row, o1, o2 in bad[:show]:
        print(f"    n={row['n']} SNF={row['D'].tolist()} HNF={row['H'][np.tril_indices(3)].tolist()} "
              f"labeling={row['labeling']}  enumlib labdegn={row['lab_degen']}  orbit={o1}  with exchange={o2}")
        print(f"      row: {row['line']}")
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("enumx", help="path to enumlib's enum.x")
    ap.add_argument("--keep", help="keep inputs and outputs in this folder")
    ap.add_argument("--show", type=int, default=3, help="example rows to print per case")
    a = ap.parse_args()
    enumx = os.path.abspath(os.path.expanduser(a.enumx))
    src = os.path.dirname(enumx)
    try:
        rev = subprocess.run(["git", "-C", src, "log", "-1", "--format=%h %cd %s", "--date=short"],
                             capture_output=True, text=True).stdout.strip()
    except OSError:
        rev = ""
    print(f"enum.x: {enumx}")
    print(f"enumlib commit: {rev or 'unknown'}")
    base = os.path.abspath(a.keep) if a.keep else tempfile.mkdtemp(prefix="labdegn_")
    print(f"{'case':22s} {'rows':>6s} {'zero':>6s} {'= orbit':>8s} {'= with ex.':>10s} "
          f"{'neither':>8s} {'Tot ok':>8s}")
    for i, (name, parent, k, sizes, mode, conc) in enumerate(CASES):
        run_case(enumx, os.path.join(base, f"case{i:02d}"), name, parent, k, sizes, mode, conc, a.show)
    print("columns: rows read; labdegn = 0; labdegn = orbit under translations x point operations;")
    print("         labdegn = orbit also under species permutations; matches neither; Totdegn = Hdegn*labdegn")
    if not a.keep:
        shutil.rmtree(base, ignore_errors=True)
    else:
        print(f"inputs and outputs kept in {base}")


if __name__ == "__main__":
    sys.exit(main())
