"""Read and write enumlib's struct_enum.out format.

The output can be used wherever enumlib's is: enumlib's makeStr.py / makestr.x
turn it into POSCAR files, and tools built on enumlib (such as pymatgen's
enumeration interface, which calls makestr.x) read it unchanged.

Each row stores a structure as superlattice + labeling, with enumlib's columns:

    #tot  HNF#  Hdegn  labdegn  Totdegn  #size  n  pg  SNF(3)  HNF(6)  L(9)  labeling

* Hdegn   number of HNFs symmetrically equivalent to this one (superlattice degeneracy)
* labdegn number of distinct labelings equivalent to this one on this superlattice
          (orbit under translations and the operations that fix the superlattice)
* Totdegn Hdegn * labdegn
* pg      number of parent space-group operations that fix the superlattice
* SNF     diagonal of the Smith normal form D = L H R
* HNF     H11 H21 H22 H31 H32 H33 (lower triangular, columns are supercell vectors
          in lattice coordinates)
* L       left transform, row by row: the site at lattice point z has address
          g = L z mod D
* labeling one digit per site, site kind first, then g in mixed-radix order:
          index = i*n + g1*d2*d3 + g2*d3 + g3

Species follow enumlib's conventions: in a "partial" listing (every species
present, no fixed composition) structures related by a permutation of the
species are identified, as in Hart & Forcade's tables; in a "full" listing and
at fixed composition they are listed separately. labdegn counts labelings
equivalent under the same group, so Totdegn summed over a cell size equals the
number of (superlattice, labeling) pairs it represents.
"""
from __future__ import annotations

import itertools

import numpy as np

from .core import _composition_vector, _label_group, _run, _site_perm, _stabilizer_ops, _sizes
from .lattice import Quotient, stabilizer, unique_superlattices
from .parent import Parent

HEADER = ("start   #tot      HNF     Hdegn   labdegn   Totdegn   #size idx    pg    SNF"
          "             HNF                 Left transform                          labeling")


def _g15(x):
    """Fortran g15.8-like formatting for a real number."""
    return f"{x:15.8f}"


def _header(P: Parent, nmin, nmax, eps, conc_ranges, include_incomplete, title):
    out = []
    out.append(' # <fortpy version="3" revision="247"></fortpy>')
    out.append(f"{title[:10]:10s}")
    out.append("bulk")
    for i in range(3):
        v = P.A[:, i]
        out.append(f"{_g15(v[0])} {_g15(v[1])} {_g15(v[2])}    # a{i + 1} parent lattice vector")
    out.append(f"{P.d:5d} # Number of points in the multilattice")
    for i, b in enumerate(P.basis):
        labs = "/".join(str(x) for x in sorted(P.allowed[i]))
        out.append(f"{_g15(b[0])} {_g15(b[1])} {_g15(b[2])}    # d{i + 1:02d} d-vector, labels: {labs}")
    out.append(f"{P.k:2d}-nary case")
    out.append(f"{nmin:4d}{nmax:4d} # Starting and ending cell sizes for search")
    out.append(f"{eps:15.8E} # Epsilon (finite precision parameter)")
    out.append("Concentration check:")
    if conc_ranges is None:
        out.append("    F")
    else:
        out.append("    T")
        out.append("Including only structures of which the concentration   of each atom is in the range:")
        for i, (lo, hi, den) in enumerate(conc_ranges):
            out.append(f"Type {i + 1}: {lo:4d}/{den:4d} -- {hi:4d}/{den:4d}")
    if include_incomplete:
        out.append("full list of labelings (including incomplete labelings) is used")
    else:
        out.append("partial list of labelings (complete labelings only) is used")
    out.append("(Non)Equivalency list:" + "".join(f"{i + 1:3d} " for i in range(P.d)))
    out.append(HEADER)
    return out


class _HNFInfo:
    """Per-superlattice data needed for the output rows."""

    def __init__(self, P: Parent, H, number):
        self.q = q = Quotient(H)
        self.H = H
        ops = _stabilizer_ops(P, H)
        pg = len(ops)
        hdegn = len(P.rots) // len(stabilizer(H, P.rots))
        n, d = q.n, P.d
        tp = q.translation_perms()
        perms = []
        for op in ops:
            p = _site_perm(P, q, op)
            for t in range(n):
                full = np.concatenate([i * n + tp[t] for i in range(d)])
                perms.append(p[full])
        G = np.unique(np.array(perms, dtype=np.int64), axis=0)
        # image of labeling lab under perm p: img[p[s]] = lab[s], i.e. img = lab[inverse(p)]
        self.inv = np.argsort(G, axis=1)
        self.hdegn = hdegn
        D, L = q.d, np.asarray(q.L)
        Hm = np.asarray(H)
        # columns that are the same for every structure on this superlattice
        self.mid = "{:3d}  {:4d}  {:2d} {:2d} {:2d}   {:2d} {:2d} {:2d} {:2d} {:2d} {:2d}   {}   ".format(
            n, pg, D[0], D[1], D[2], int(Hm[0, 0]), int(Hm[1, 0]), int(Hm[1, 1]), int(Hm[2, 0]),
            int(Hm[2, 1]), int(Hm[2, 2]), " ".join(f"{int(x):4d}" for x in L.ravel()))
        self.number = number

    def labdegn(self, labs, LG):
        """Orbit sizes of a batch of labelings (rows of labs) under this superlattice's
        translations and fixing operations, and the label permutations LG:
        |group| / |stabilizer|, found by applying every site permutation of the group.
        The writer takes its degeneracies from the generator; this independent
        calculation is what tests/test_degeneracy.py compares them with."""
        stab = np.zeros(len(labs), dtype=np.int64)
        for e in LG:
            el = e[labs]
            for inv in self.inv:
                stab += (el[:, inv] == labs).all(axis=1)
        return (len(self.inv) * len(LG)) // stab


def _compositions_in_ranges(P: Parent, N, ranges, include_incomplete):
    """Atom counts c (summing to N) with lo_i/den_i <= c_i/N <= hi_i/den_i for every species."""
    out = []
    for cut in itertools.combinations_with_replacement(range(N + 1), P.k - 1):
        c = np.diff((0,) + cut + (N,))
        if not include_incomplete and c.min() == 0:
            continue
        if all(lo * N <= ci * den <= hi * N for ci, (lo, hi, den) in zip(c, ranges)):
            out.append(tuple(int(x) for x in c))
    return out


def write_struct_enum(file, parent: Parent, sizes, composition=None, include_incomplete=False,
                      title="necklace-enum", eps=1e-6, concentration_ranges=None):
    """Enumerate and write the structures in enumlib's struct_enum.out format.

    Parameters
    ----------
    file : path or text file object
    parent, sizes, composition, include_incomplete
        As for :func:`necklace_enum.enumerate_structures`. A composition given as
        atom counts needs a single size.
    concentration_ranges : list of (lo, hi, den), optional
        enumlib-style ranges, one per species: lo/den <= fraction <= hi/den.
        Every composition in range is enumerated for each size.
    title : str
        Up to 10 characters, written on the title line.

    Label exchange follows enumlib: it is applied in a partial listing (every
    species present) without concentration restrictions, and not otherwise.
    Returns the number of structures written.
    """
    if parent.k > 10:
        raise ValueError("the struct_enum.out format stores one digit per site (at most 10 species)")
    sizes = sorted(_sizes(sizes))
    if composition is not None and concentration_ranges is not None:
        raise ValueError("give either composition or concentration_ranges, not both")
    if composition is not None:
        if len(sizes) != 1:
            raise ValueError("a fixed composition (atom counts) needs a single size; use concentration_ranges")
        N = sizes[0] * parent.d
        c = _composition_vector(parent, composition, N)
        if c is None:
            raise ValueError(f"composition does not add up to {N} sites")
        concentration_ranges = [(int(x), int(x), N) for x in c]
    if concentration_ranges is not None and len(concentration_ranges) != parent.k:
        raise ValueError("need one concentration range per species")
    label_exchange = concentration_ranges is None and not include_incomplete
    LG = _label_group(parent, label_exchange)
    own = isinstance(file, (str, bytes)) or hasattr(file, "__fspath__")
    fh = open(file, "w") if own else file
    try:
        for line in _header(parent, sizes[0], sizes[-1], eps, concentration_ranges, include_incomplete, title):
            fh.write(line + "\n")
        tot = 0
        hnf_offset = 0
        for n in sizes:
            reps = unique_superlattices(n, parent.rots)
            hnf_number = {tuple(np.asarray(H).ravel()): hnf_offset + i + 1 for i, H in enumerate(reps)}
            hnf_offset += len(reps)
            if concentration_ranges is None:
                runs = [None]
            else:
                runs = _compositions_in_ranges(parent, n * parent.d, concentration_ranges, include_incomplete)
            info = {}
            size_count = 0
            for comp in runs:
                for Hs, labs, degs, _ in _run(parent, n, "batch", label_exchange, comp, include_incomplete, None):
                    strings = (labs.astype(np.uint8) + 48).view(f"S{labs.shape[1]}").ravel()
                    for H in Hs:
                        key = tuple(np.asarray(H).ravel())
                        h = info.get(key)
                        if h is None:
                            h = info[key] = _HNFInfo(parent, np.asarray(H), hnf_number[key])
                        lds = degs                     # from the generator's symmetry comparisons
                        rows = []
                        for ld, lab in zip(lds.tolist(), strings):
                            tot += 1
                            size_count += 1
                            rows.append(f"{tot:11d} {h.number:9d} {h.hdegn:7d} {ld:8d} {ld * h.hdegn:8d} "
                                        f"{size_count:11d} {h.mid}{lab.decode()}\n")
                        fh.write("".join(rows))
        return tot
    finally:
        if own:
            fh.close()


def read_struct_enum_in(file):
    """Read an enumlib input file (struct_enum.in).

    Returns a dict with 'parent' (a Parent; species are labelled 0..k-1),
    'sizes', 'eps', 'include_incomplete' (True for a "full" listing),
    'concentration_ranges' (None or a list of (lo, hi, den)) and 'title'.
    """
    own = isinstance(file, (str, bytes)) or hasattr(file, "__fspath__")
    fh = open(file) if own else file
    try:
        raw = [l.split("#")[0].strip() for l in fh]
    finally:
        if own:
            fh.close()
    lines = [l for l in raw if l]
    title = lines[0]
    if not lines[1].lower().startswith("b"):
        raise ValueError("only bulk (3D) parents are supported")
    lat = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    k = int(lines[5].split("-")[0].split()[0])
    nd = int(lines[6].split()[0])
    basis, allowed = [], []
    for i in range(nd):
        f = lines[7 + i].split()
        basis.append([float(x) for x in f[:3]])
        allowed.append([int(x) for x in f[3].split("/")] if len(f) > 3 else list(range(k)))
    j = 7 + nd
    if lines[j].lower().startswith("e"):            # optional site-equivalency line
        j += 2
    nmin, nmax = (int(x) for x in lines[j].split()[:2])
    eps = float(lines[j + 1].split()[0].replace("D", "E").replace("d", "e"))
    mode = lines[j + 2].lower()
    if not (mode.startswith("full") or mode.startswith("part")):
        raise ValueError(f"expected 'full' or 'part', found: {lines[j + 2]!r}")
    ranges = None
    rest = lines[j + 3:]
    if rest:                                        # one "min max denominator" line per species
        if len(rest) < k:
            raise ValueError(f"concentration ranges need {k} lines (one per species), found {len(rest)}")
        ranges = [tuple(int(x) for x in rest[i].split()[:3]) for i in range(k)]
        for lo, hi, den in ranges:
            if not 0 <= lo <= hi <= den:
                raise ValueError(f"bad concentration range {lo} {hi} {den}")
    parent = Parent(lat, basis, species=k, allowed=allowed)
    return {"parent": parent, "sizes": list(range(nmin, nmax + 1)), "eps": eps,
            "include_incomplete": mode.startswith("full"), "concentration_ranges": ranges, "title": title}


def run_enumlib_input(infile="struct_enum.in", outfile="struct_enum.out"):
    """Read an enumlib input file and write the matching struct_enum.out, as enum.x does."""
    spec = read_struct_enum_in(infile)
    return write_struct_enum(outfile, spec["parent"], spec["sizes"], include_incomplete=spec["include_incomplete"],
                             title=spec["title"], eps=spec["eps"], concentration_ranges=spec["concentration_ranges"])


def read_struct_enum(file):
    """Read a struct_enum.out file (from enumlib or necklace-enum).

    Returns (header, rows): header is a dict with 'title', 'lattice' (rows are the
    parent vectors), 'basis', 'k'; rows is a list of dicts with keys 'n', 'hnf'
    (3x3), 'snf' (3,), 'L' (3x3), 'labeling' (str), 'hnf_degen', 'lab_degen',
    'tot_degen', 'pg'.
    """
    own = isinstance(file, (str, bytes)) or hasattr(file, "__fspath__")
    fh = open(file) if own else file
    try:
        lines = [l.rstrip("\n") for l in fh]
    finally:
        if own:
            fh.close()
    body = [l for l in lines if l.strip() and not l.lstrip().startswith("#")]
    header = {"title": body[0].strip()}
    header["lattice"] = np.array([[float(x) for x in body[i].split()[:3]] for i in (2, 3, 4)])
    nd = int(body[5].split()[0])
    header["basis"] = np.array([[float(x) for x in body[6 + i].split()[:3]] for i in range(nd)])
    header["k"] = int(body[6 + nd].split("-")[0])
    start = next(i for i, l in enumerate(lines) if l.startswith("start"))
    rows = []
    for l in lines[start + 1:]:
        f = l.split()
        if len(f) < 27:
            continue
        v = [int(x) for x in f[:26]]
        H = np.array([[v[11], 0, 0], [v[12], v[13], 0], [v[14], v[15], v[16]]])
        rows.append({"n": v[6], "hnf_degen": v[2], "lab_degen": v[3], "tot_degen": v[4], "pg": v[7],
                     "snf": np.array(v[8:11]), "hnf": H, "L": np.array(v[17:26]).reshape(3, 3),
                     "labeling": f[26]})
    return header, rows
