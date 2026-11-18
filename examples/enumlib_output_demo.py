#!/usr/bin/env python3
"""Demonstration of necklace-enum's enumlib-compatible input and output.

Run from the package folder:

    python examples/enumlib_output_demo.py

Everything is written to ./demo_output/. Needs the package (pip install .);
ASE is optional (pip install ase) for the last step.
"""
import collections
import os
import sys

try:
    import necklace_enum as te
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
    import necklace_enum as te

OUT = "demo_output"
os.makedirs(OUT, exist_ok=True)


def rule(title):
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78)


# --------------------------------------------------------------------------------------
rule("1. Write a struct_enum.out file: binary fcc, cells of 1-6 atoms")
P = te.fcc(2)
path = os.path.join(OUT, "fcc_binary_struct_enum.out")
n = te.write_struct_enum(path, P, range(1, 7))
print(f"wrote {path}: {n} structures\n")
lines = open(path).read().splitlines()
start = next(i for i, l in enumerate(lines) if l.startswith("start"))
print("Header (same layout as enumlib's):")
for l in lines[:start + 1]:
    print("   " + l)
print("\nFirst 8 structures:")
for l in lines[start + 1:start + 9]:
    print("   " + l)
print("   ...")

# --------------------------------------------------------------------------------------
rule("2. What one row means")
header, rows = te.read_struct_enum(path)
r = next(r for r in rows if r["n"] == 4 and r["labeling"].count("1") == 2)
print(f"""Structure #{rows.index(r) + 1}, a 4-atom cell with labeling {r['labeling']}:
   n (atoms per cell)           {r['n']}
   superlattice (HNF)           {r['hnf'].tolist()}
   Smith normal form diagonal   {r['snf'].tolist()}
   left transform L             {r['L'].tolist()}
   labeling                     {r['labeling']}   (one digit per site: 0 = A, 1 = B)
   pg                           {r['pg']}   operations of the parent lattice keep this supercell fixed
   Hdegn                        {r['hnf_degen']}   symmetry-equivalent supercells
   labdegn                      {r['lab_degen']}   equivalent labelings on this supercell
   Totdegn                      {r['tot_degen']}   = Hdegn x labdegn configurations in total
The atom at lattice point z sits at group address g = L z mod D, and its type is the
digit of the labeling at position  g1*d2*d3 + g2*d3 + g3  (enumlib's rule).""")

# --------------------------------------------------------------------------------------
rule("3. Degeneracies add up")
by_n = collections.defaultdict(lambda: [0, 0])
for r in rows:
    by_n[r["n"]][0] += 1
    by_n[r["n"]][1] += r["tot_degen"]
print("For each cell size, Totdegn summed over the structures counts every (supercell,")
print("labeling) pair they represent, as a cluster-expansion or partition-function sum needs.\n")
print(f"   {'n':>3} {'structures':>11} {'sum of Totdegn':>15}")
for k_, (cnt, tot) in sorted(by_n.items()):
    print(f"   {k_:>3} {cnt:>11} {tot:>15}")

# --------------------------------------------------------------------------------------
rule("4. Drop-in for enumlib's enum.x: read a struct_enum.in, write struct_enum.out")
inp = """fcc
bulk
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5
 2 -nary case
 1  # Number of points in the multilattice
 0.0 0.0 0.0  0/1
 1 8   # Starting and ending cell sizes for search
0.10000000E-06
full list of labelings
# Concentration restrictions: fraction of A between 1/2 and 1, of B between 0 and 1
1 2 2
0 2 2
"""
d = os.path.join(OUT, "dropin")
os.makedirs(d, exist_ok=True)
open(os.path.join(d, "struct_enum.in"), "w").write(inp)
print("struct_enum.in (enumlib's own example input):")
for l in inp.strip().splitlines():
    print("   " + l)
n = te.run_enumlib_input(os.path.join(d, "struct_enum.in"), os.path.join(d, "struct_enum.out"))
_, rows2 = te.read_struct_enum(os.path.join(d, "struct_enum.out"))
print(f"\nwrote {d}/struct_enum.out: {n} structures (enumlib's enum.x also writes 376)")
c = collections.Counter(r["n"] for r in rows2)
print("   per cell size:", dict(sorted(c.items())))
print("   on the command line:  cd", d, "&& necklace-enum enumlib")

# --------------------------------------------------------------------------------------
rule("5. Other parents: hcp with a restricted site")
H = te.hcp(["Li", "Vac", "Co"], allowed=[{"Li", "Vac"}, {"Co", "Vac"}])
path3 = os.path.join(OUT, "hcp_restricted_struct_enum.out")
n = te.write_struct_enum(path3, H, range(1, 4), include_incomplete=True)
lines = open(path3).read().splitlines()
print(f"wrote {path3}: {n} structures")
print("d-vector lines record which species may sit on each site:")
for l in lines:
    if "d-vector" in l:
        print("   " + l)

# --------------------------------------------------------------------------------------
rule("6. From a row back to a crystal")
try:
    import ase  # noqa: F401
    s = next(s for s in te.enumerate_structures(te.fcc(["Cu", "Au"]), 4, composition={"Cu": 2, "Au": 2}))
    atoms = s.to_ase()
    print(f"A {len(atoms)}-atom Cu2Au2 structure built by necklace-enum:")
    print("   cell (rows):", atoms.cell[:].round(3).tolist())
    for sym, pos in zip(atoms.get_chemical_symbols(), atoms.positions.round(3).tolist()):
        print(f"   {sym:2s} {pos}")
    vasp = os.path.join(OUT, "POSCAR_Cu2Au2")
    atoms.write(vasp, format="vasp")
    print(f"written as {vasp}")
except ImportError:
    print("(install ASE to build crystals: pip install ase)")
print("enumlib's makeStr.py or makestr.x can also turn any struct_enum.out written here into POSCAR files.")

print(f"\nAll output is in ./{OUT}/")
