from pathlib import Path

from structphylogeny.io import choose_best_chain


def test_choose_best_chain_extracts_sequence():
    record = choose_best_chain(
        Path("data/mouse_lipocalin_structures/APOD.cif")
    )
    assert record.sample == "APOD"
    assert record.chain_id == "A0"
    assert record.residue_count == len(record.sequence)
    assert record.residue_count > 100
    assert record.sequence.startswith("QNFHLGKCPS")
