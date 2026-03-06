from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1


STANDARD_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "MSE",
}


@dataclass(frozen=True)
class StructureRecord:
    sample: str
    path: Path
    chain_id: str
    sequence: str
    residue_count: int


def iter_structure_files(structures_dir: Path) -> list[Path]:
    suffixes = {".pdb", ".ent", ".cif", ".mmcif"}
    files = [
        path
        for path in sorted(structures_dir.iterdir())
        if path.is_file() and any(path.name.endswith(ext) for ext in suffixes)
    ]
    if not files:
        raise FileNotFoundError(f"No supported structure files found in {structures_dir}")
    return files


def _load_structure(path: Path):
    if path.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def _residue_to_aa(resname: str) -> str:
    return seq1(resname, custom_map={"MSE": "M", "SEC": "U", "PYL": "O"})


def _chain_sequence(chain) -> str:
    letters: list[str] = []
    for residue in chain.get_residues():
        if residue.id[0] != " ":
            continue
        if residue.resname not in STANDARD_RESIDUES:
            continue
        aa = _residue_to_aa(residue.resname)
        if aa == "X":
            continue
        letters.append(aa)
    return "".join(letters)


def choose_best_chain(path: Path) -> StructureRecord:
    structure = _load_structure(path)
    chains = []
    for chain in structure.get_chains():
        sequence = _chain_sequence(chain)
        if not sequence:
            continue
        chains.append((chain.id, sequence))
    if not chains:
        raise ValueError(f"No protein chain with standard residues found in {path}")
    chain_id, sequence = max(chains, key=lambda item: (len(item[1]), item[0]))
    sample = path.name
    for suffix in (".cif", ".mmcif", ".pdb", ".ent"):
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
            break
    return StructureRecord(
        sample=sample,
        path=path,
        chain_id=chain_id,
        sequence=sequence,
        residue_count=len(sequence),
    )


def load_records(structures_dir: Path) -> list[StructureRecord]:
    return [choose_best_chain(path) for path in iter_structure_files(structures_dir)]


def write_fasta(records: Iterable[StructureRecord], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.sample}\n{record.sequence}\n")
