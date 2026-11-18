# necklace-enum

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24330096.svg)](https://doi.org/10.5281/zenodo.24330096)

Enumerates derivative superstructures — the symmetrically distinct ways of decorating a parent lattice with
several kinds of atoms — by generating them directly as canonical words, without filtering duplicates.

Project page: https://r-h-taylor.github.io/necklace-enum/ · Archived at Zenodo: https://doi.org/10.5281/zenodo.24330096

This is the reference implementation for

> R. Taylor, *Improving the Hart-Forcade enumeration algorithm: derivative superstructures as Lyndon words* (2026).

It keeps the superlattice and Smith-normal-form construction of Hart & Forcade
([Phys. Rev. B 77, 224115 (2008)](https://doi.org/10.1103/PhysRevB.77.224115)) and replaces their duplicate-removal
step. Reading the sites column by column along the largest cyclic factor of the quotient group turns each labeling
into a word. Translations act on the word as rotations combined with a single letterwise permutation, so the
canonical structures are Lyndon words, which are generated directly with an extended necklace recursion.

## Install

```bash
pip install .            # numpy + numba
pip install ".[ase]"     # to export ASE Atoms
pip install ".[test]"    # to run the tests
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

**Conventions.**
- By default, species are distinct: A₃B and AB₃ are different structures. With `label_exchange=True`, structures
  that differ only by a permutation of species are identified, as in Hart & Forcade's tables.
- Only structures that contain every species are listed unless `include_incomplete=True`.
- Structures that are periodic in a smaller cell are excluded; they appear at their own, smaller size.

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
tabulated by Hart & Forcade, and check that each exported crystal is primitive in its own cell and that no two
exported crystals are symmetry-equivalent. `tests/hart_forcade.py` is a compiled reference implementation of the 2008 algorithm, used to check counts and to compare speed.

## Performance and limits

- The compiled generator labels binary fcc structures at tens of nanoseconds per structure. It lists all
  4,592,213,143 binary fcc structures with 32 sites per cell in about six minutes on one core (`count_structures`).
- `enumerate_structures` builds a Python object per structure, so listing is far slower than counting. Use counts,
  fixed compositions or streaming for very large runs.
- The per-superlattice setup is plain Python and takes a few tenths of a second per size. For dilute fixed-composition
  runs it is larger than the labeling time.
- Letters hold `d * d1 * d2` sites. For multilattices with many sites and species the letter tables can grow large,
  in which case a `MemoryError` names the offending Smith normal form.
- The first call compiles the kernels with Numba (a few seconds); later calls use the on-disk cache.

## License

MIT. If you use this code, please cite the paper above and Hart & Forcade (2008).
