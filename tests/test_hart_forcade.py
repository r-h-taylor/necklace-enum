"""The Hart-Forcade reference implementation agrees with necklace-enum."""
import pytest

import necklace_enum as te

import hart_forcade as hf


@pytest.mark.parametrize("k,nmax", [(2, 12), (3, 8), (4, 7)])
def test_reference_agrees(k, nmax):
    P = te.fcc(k)
    for n in range(2, nmax + 1):
        assert hf.count_structures(P, n) == te.count_structures(P, n, label_exchange=True)[n]
