"""enumerate_labelings (array listing) gives exactly the structures of enumerate_structures."""
import numpy as np
import pytest

import necklace_enum as te
from necklace_enum import core

CASES = [
    (lambda: te.fcc(2), range(1, 11), None, True),
    (lambda: te.fcc(2), range(1, 11), None, False),
    (lambda: te.fcc(3), range(1, 7), None, False),
    (lambda: te.fcc(3), 6, (2, 2, 2), False),
    (lambda: te.fcc(4), range(1, 6), None, False),
    (lambda: te.bcc(2), range(1, 9), None, False),
    (lambda: te.sc(3), range(1, 5), None, True),
    (lambda: te.hcp(2), range(1, 5), None, False),
    (lambda: te.hcp(3), range(1, 4), None, False),
    (lambda: te.fcc(2), 12, (9, 3), False),
]


def as_set(pairs):
    out = set()
    for H, lab in pairs:
        out.add((np.asarray(H).tobytes(), tuple(int(x) for x in lab)))
    return out


@pytest.mark.parametrize("make,sizes,comp,lx", CASES)
@pytest.mark.parametrize("batch", [3, 1 << 16])
def test_labelings_match_structures(make, sizes, comp, lx, batch, monkeypatch):
    monkeypatch.setattr(core, "LIST_BATCH", batch)
    P = make()
    ref = as_set((s.hnf, s.labels) for s in te.enumerate_structures(P, sizes, composition=comp, label_exchange=lx))
    rows = []
    for H, labs in te.enumerate_labelings(P, sizes, composition=comp, label_exchange=lx):
        assert labs.ndim == 2 and not labs.flags.writeable
        rows.extend((H, r) for r in labs)
    got = as_set(rows)
    assert len(rows) == len(got) == len(ref)            # no duplicates, nothing missing
    assert got == ref
    cnt = te.count_structures(P, sizes, composition=comp, label_exchange=lx)
    assert len(rows) == sum(cnt.values())


def test_rebuild_structure():
    P = te.fcc(["Cu", "Au"])
    H, labs = next(iter(te.enumerate_labelings(P, 4)))
    s = te.Structure(P, H, tuple(int(x) for x in labs[0]))
    assert len(s.sites()[0]) == 4
