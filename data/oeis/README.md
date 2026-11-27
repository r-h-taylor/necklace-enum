# OEIS submission drafts: binary fcc derivative structures

Two companion sequences, n = 1..100, computed with `count_structures_burnside` (`scripts/oeis_counts.py`):

| File | Sequence |
|---|---|
| `b_fcc_binary_label_exchange.txt` | up to interchange of the two species (Hart & Forcade's convention) |
| `b_fcc_binary.txt` | the two species distinguished |

**Checks.** Every term with n <= 32 equals the count from direct generation (`count_structures`), which is itself
tested structure by structure against brute force. Terms n = 2..23 of the first sequence equal Table V of Hart &
Forcade (2008). The superlattice counts used by both codes equal A159842 for n = 1..60. The fcc and bcc lattices
give identical counts (the same holds for A159842).

Before submitting, search oeis.org for the first dozen terms of each sequence to confirm they are not already in
the database. Draft entries follow in OEIS internal format; replace `A??????` by the assigned numbers.

---

## Entry 1: up to interchange of the species

```
%N Number of binary derivative superstructures of the face-centered cubic lattice with n sites per primitive cell, up to the symmetry of the lattice and interchange of the two atom types.
%D 0,2,3,12,14,50,52,229,252,685,682,3875,2624,9628,16584,49764,42135,212612,174104,867893,1120708,2628180,3042732,19335101,15592901,45995141,82823587,255215928,224864930,1311471651,941857818,4592213143,5925530526,14092259460
%O 1,2
%C A derivative superstructure is a decoration of the lattice by two kinds of atom, A and B, that repeats with a supercell of n primitive cells and has no smaller period. Both kinds must occur. Two decorations are counted once if a symmetry of the lattice (rotation, reflection or translation), possibly combined with interchange of A and B, maps one onto the other.
%C Decorations that repeat with a smaller cell are counted at their own smaller n, so a(1) = 0.
%C The same sequence counts binary derivative superstructures of the body-centered cubic lattice.
%C For each of the A159842(n) symmetry-distinct superlattices S of index n, the structures are the orbits of labelings of the quotient group G = L/S (of order n) that are fixed by no nonzero translation, under translations, the point operations fixing S, and interchange of the labels. Burnside's lemma counts these orbits; Moebius inversion over the subgroups of G removes the labelings with a smaller period. Only subgroups whose Sylow subgroups are elementary abelian contribute.
%H Richard H. Taylor, <a href="/A??????/b??????.txt">Table of n, a(n) for n = 1..100</a>
%H Gus L. W. Hart and Rodney W. Forcade, <a href="https://doi.org/10.1103/PhysRevB.77.224115">Algorithm for generating derivative structures</a>, Phys. Rev. B 77, 224115 (2008); <a href="https://arxiv.org/abs/0804.3544">arXiv:0804.3544</a>. Table V gives a(2)..a(23).
%H Richard H. Taylor, <a href="https://github.com/r-h-taylor/necklace-enum">necklace-enum</a>, derivative structures as Lyndon words (Python), 2026; <a href="https://doi.org/10.5281/zenodo.24330095">doi:10.5281/zenodo.24330095</a>.
%e a(2) = 2: the L1_0 (CuAu) structure, alternating (001) planes of A and B, and the L1_1 (CuPt) structure, alternating (111) planes.
%e a(3) = 3: one A2B structure on each of the A159842(3) = 3 superlattices of index 3.
%o (Python) # pip install git+https://github.com/r-h-taylor/necklace-enum
%o import necklace_enum as te
%o def a(n): return te.count_structures_burnside(te.fcc(2), n, label_exchange=True)[n]
%Y Cf. A001001, A159842, A?????? (species distinguished).
%K nonn
%A Richard H. Taylor, Nov 2026
```

## Entry 2: species distinguished

```
%N Number of binary derivative superstructures of the face-centered cubic lattice with n sites per primitive cell, up to the symmetry of the lattice, with the two atom types distinguished.
%D 0,2,6,19,28,80,104,390,504,1211,1364,7140,5248,18270,33168,95837,84270,417268,348208,1715551,2241416,5226228,6085464,38517429,31185802,91834053,165647174,509897374,449729860,2621617654,1883715636,9181330073,11851061052
%O 1,2
%C As A??????, but a decoration and the one obtained by interchanging A and B are counted separately unless a symmetry of the lattice maps one onto the other. For example, A3B and AB3 structures are always distinct.
%C The same sequence counts binary derivative superstructures of the body-centered cubic lattice.
%C A??????(n) <= a(n) <= 2*A??????(n).
%H Richard H. Taylor, <a href="/A??????/b??????.txt">Table of n, a(n) for n = 1..100</a>
%H Gus L. W. Hart and Rodney W. Forcade, <a href="https://doi.org/10.1103/PhysRevB.77.224115">Algorithm for generating derivative structures</a>, Phys. Rev. B 77, 224115 (2008).
%H Richard H. Taylor, <a href="https://github.com/r-h-taylor/necklace-enum">necklace-enum</a>, derivative structures as Lyndon words (Python), 2026.
%e a(2) = 2: L1_0 and L1_1 each have one A and one B per cell, and a translation interchanges them.
%e a(3) = 6: each of the three A2B structures and its AB2 counterpart.
%o (Python) # pip install git+https://github.com/r-h-taylor/necklace-enum
%o import necklace_enum as te
%o def a(n): return te.count_structures_burnside(te.fcc(2), n)[n]
%Y Cf. A001001, A159842, A?????? (up to interchange of species).
%K nonn
%A Richard H. Taylor, Nov 2026
```
