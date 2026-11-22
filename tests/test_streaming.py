"""Listing streams through a fixed buffer: with tiny batches it must give exactly the
structures of a single large batch, in the same order, and agree with the counter."""
import pytest

import necklace_enum as te
import necklace_enum.core as core

CASES = [(te.fcc(2), 8, None, True), (te.fcc(3), 6, None, False), (te.hcp(2), 4, None, False),
         (te.fcc(2), 12, (9, 3), False), (te.bcc(2), 7, None, False)]


def listing(P, n, comp, le):
    return [(s.hnf.tobytes(), s.labels) for s in te.enumerate_structures(P, n, composition=comp, label_exchange=le)]


@pytest.mark.parametrize("P,n,comp,le", CASES)
@pytest.mark.parametrize("batch", [1, 2, 5, 64])
def test_small_batches_match(P, n, comp, le, batch, monkeypatch):
    ref = listing(P, n, comp, le)
    monkeypatch.setattr(core, "LIST_BATCH", batch)
    got = listing(P, n, comp, le)
    assert got == ref
    assert len(got) == te.count_structures(P, n, composition=comp, label_exchange=le)[n]
