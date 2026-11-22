# necklace-enum

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24330096.svg)](https://doi.org/10.5281/zenodo.24330096)

Enumerates derivative superstructures — the symmetrically distinct ways of decorating a parent lattice with
several kinds of atoms — by generating them directly as canonical words, without filtering duplicates.

Project page: https://r-h-taylor.github.io/necklace-enum/ · Archived at Zenodo: https://doi.org/10.5281/zenodo.24330096

This is the reference implementation for

> R. H. Taylor, *Toroidal necklaces for low-memory derivative structure enumeration* (2026).

It keeps the superlattice and Smith-normal-form construction of Hart & Forcade
([Phys. Rev. B 77, 224115 (2008)](https://doi.org/10.1103/PhysRevB.77.224115)) and replaces their duplicate-removal
step. Reading the sites column by column along the largest cyclic factor of the quotient group turns each labeling
into a word. Translations act on the word as rotations combined with a single letterwise permutation, so the
canonical structures are Lyndon words, which are generated directly with an extended necklace recursion.

## Install

```bash
pip install necklace-enum            # numpy + numba
pip install "necklace-enum[ase]"     # to export ASE Atoms
```

From a clone of this repository (needed for the tests and the benchmark scripts):

```bash
pip install ".[test]"
pytest
```

## Use

```python
import necklace_enum as te

P = te.fcc(["Cu", "Au"])
te.count_structures(P, range(1, 13))                 # {n: count} without listing them

for s in te.enumerate_structures(P, 8, composition={"Cu": 6, "Au": 2}):
    atoms = s.to_ase()                               # or s.to_pymatgen(), or s.hnf / s.labels

hcp = te.hcp(["Li", "Vac", "Co"], allowed=[{"Li", "Vac"}, {"Co", "Vac"}])   # site restrictions
custom = te.Parent(lattice, basis, species=["A", "B", "C"])                 # any parent
```

Command line:

```bash
necklace-enum count fcc Cu,Au 1-12
necklace-enum count fcc 2 2-23 --label-exchange          # Hart & Forcade, Table V
necklace-enum list hcp Mg,Cd 4 --composition 5,3 --format extxyz -o out.xyz
```

**Counting without listing.** `count_structures_burnside` gives the same counts as `count_structures`, with the
same arguments, by Burnside's lemma: it sums fixed points over the symmetry operations of each superlattice instead
of generating structures. Its cost grows with the number of superlattices, not structures, so it reaches far larger
cells (binary fcc at n = 32 in about a second, n = 96 in under half a minute):

```python
te.count_structures_burnside(te.fcc(2), range(1, 61), label_exchange=True)
```

```bash
necklace-enum count fcc 2 1-60 --label-exchange --burnside
python scripts/oeis_counts.py --nmax 100 -o b_fcc_binary.txt     # OEIS b-file
```

**Conventions.**
- By default, species are distinct: A₃B and AB₃ are different structures. With `label_exchange=True`, structures
  that differ only by a permutation of species are identified, as in Hart & Forcade's tables.
- Only structures that contain every species are listed unless `include_incomplete=True`.
- Structures that are periodic in a smaller cell are excluded; they appear at their own, smaller size.

## enumlib compatibility

necklace-enum reads enumlib's input file and writes its `struct_enum.out` format, so it can take the place of
`enum.x` in existing workflows:

```bash
necklace-enum enumlib struct_enum.in                    # writes struct_enum.out, as enum.x does
necklace-enum list fcc Cu,Au 1-12 --format enum -o struct_enum.out
```

```python
te.run_enumlib_input("struct_enum.in", "struct_enum.out")
te.write_struct_enum("struct_enum.out", P, range(1, 13))
header, rows = te.read_struct_enum("struct_enum.out")    # also reads enumlib's own files
```

enumlib's `makeStr.py` and `makestr.x` (which pymatgen calls) build POSCAR files from this output unchanged.
The output follows enumlib's conventions:
- A "partial" listing without concentration restrictions identifies structures related by swapping species, as in
  Hart & Forcade's tables. "Full" listings and concentration-restricted runs list them separately.
- The degeneracy columns are exact: `Totdegn` summed over a cell size equals the number of (superlattice, labeling)
  pairs the structures stand for. enumlib writes `labdegn = 0` in its fixed-concentration mode.

`scripts/compare_with_enumlib.py` runs both programs on the same input and compares the results.

## Verification

```bash
python check_necklace_enum.py           # stand-alone check incl. timing against Hart & Forcade (under a minute)
python check_necklace_enum.py --full    # more cases and larger cells (a few minutes)
pytest                   # brute force, published counts, export checks (about 15 s)
pytest --slow            # also binary fcc up to n = 23
python scripts/reproduce_paper.py [--full]
```

The tests compare the generated structures, one by one, with an independent brute-force code for fcc, bcc,
simple-cubic and hcp parents (including site restrictions and fixed compositions). They also reproduce every count
tabulated by Hart & Forcade, compare the Burnside counter with the generator over parents, species counts,
compositions and site restrictions, and check that each exported crystal is primitive in its own cell and that no two
exported crystals are symmetry-equivalent. The enumlib-format tests rebuild every structure from the written file with enumlib's indexing rule and check the degeneracy columns against a direct count. `tests/hart_forcade.py` is a compiled reference implementation of the 2008 algorithm, used to check counts and to compare speed.

## Performance and limits

- The compiled generator labels binary fcc structures at tens of nanoseconds per structure. `count_structures`
  counts all 4,592,213,143 binary fcc structures with 32 sites per cell (up to label exchange) in about 230 s on
  one core of an Apple M2.
- Memory does not grow with the number of structures. The generator's state is polynomial in the cell size, and
  listing streams structures from the compiled kernel in batches of 2^16 words. Counting stays within 10 MB of the
  85 MB interpreter baseline up to 9.2 billion structures; listing stays below 0.22 GB
  (`scripts/memory_compare.py`).
- `enumerate_structures` builds a Python object per structure, so listing is far slower than counting. Use counts
  or fixed compositions for very large runs.
- The per-superlattice setup is plain Python and takes a few tenths of a second per size. For dilute fixed-composition
  runs it is larger than the labeling time.
- Letters hold `d * d1 * d2` sites. For multilattices with many sites and species the letter tables can grow large,
  in which case a `MemoryError` names the offending Smith normal form.
- The first call compiles the kernels with Numba (a few seconds); later calls use the on-disk cache.

## Changes

- **0.2.0** (21 Nov 2026): streaming listing with constant memory (resumable kernel); Burnside counter
  (`burnside.py`); benchmark and comparison scripts against enumlib, Enumlib.jl and dsenum; structure-by-structure
  check against enumlib (`scripts/compare_structures_with_enumlib.py`); binary fcc counts to 43 sites
  (`scripts/oeis_counts.py`); current paper in `docs/`.
- **0.1.0** (18 Nov 2026): first release.

## License

MIT. If you use this code, please cite the paper above and Hart & Forcade (2008).
