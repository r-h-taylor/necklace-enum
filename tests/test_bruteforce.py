"""Structure-by-structure comparison with the brute-force reference."""
from collections import defaultdict

import pytest

import necklace_enum as te
from necklace_enum.lattice import unique_superlattices

from brute import brute, canon, check_superlattices, label_group, to_cell_labeling

PARENTS = {
    "fcc2": lambda: te.fcc(2), "fcc3": lambda: te.fcc(3), "bcc2": lambda: te.bcc(2),
    "sc3": lambda: te.sc(3), "hcp2": lambda: te.hcp(2), "hcp3": lambda: te.hcp(3),
    "hcp3r": lambda: te.hcp(3, allowed=[{0, 1}, {1, 2}]),
}
CASES = (
    [("fcc2", n, le, None) for n in range(1, 8) for le in (False, True)]
    + [("fcc3", n, le, None) for n in range(1, 5) for le in (False, True)]
    + [("bcc2", n, le, None) for n in range(1, 7) for le in (False, True)]
    + [("sc3", n, le, None) for n in range(1, 4) for le in (False, True)]
    + [("hcp2", n, le, None) for n in range(1, 5) for le in (False, True)]
    + [("hcp3", n, le, None) for n in range(1, 3) for le in (False, True)]
    + [("hcp3r", n, le, None) for n in range(1, 4) for le in (False, True)]
    + [("fcc2", 8, False, (5, 3)), ("fcc3", 6, False, (2, 2, 2)), ("hcp2", 4, False, (5, 3)),
       ("hcp3", 3, False, (2, 2, 2)), ("hcp3r", 4, False, (3, 2, 3)), ("hcp3r", 4, True, (3, 2, 3))]
)


@pytest.mark.parametrize("name,n,le,conc", CASES)
def test_against_brute_force(name, n, le, conc):
    P = PARENTS[name]()
    reps = unique_superlattices(n, P.rots)
    check_superlattices(P, n, reps)
    ref = brute(P, n, reps, le, conc)
    LG = label_group(P, le)
    got = defaultdict(list)
    for s in te.enumerate_structures(P, n, composition=conc, label_exchange=le):
        cell, group, _ = ref[tuple(s.hnf.ravel())]
        got[tuple(s.hnf.ravel())].append(canon(to_cell_labeling(P, cell, s.hnf, s.labels), group, LG))
    for h, (_, _, expect) in ref.items():
        mine = got.get(h, [])
        assert len(mine) == len(set(mine)), f"duplicate structures for H={h}"
        assert set(mine) == expect, f"H={h}: {len(set(mine))} structures, brute force {len(expect)}"
    total = sum(len(e) for _, _, e in ref.values())
    assert te.count_structures(P, n, composition=conc, label_exchange=le)[n] == total
