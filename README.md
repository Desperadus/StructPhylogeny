# StructPhylogeny

`StructPhylogeny` is a Snakemake pipeline for building parallel structure-based and sequence-based phylogenies from a folder of protein structures.

It is designed around the workflow direction of Lakshmi et al. (2015), but implemented with modern tooling:

- `GTalign` for all-vs-all 3D structural alignment and comparison
- an SDM-like structural dissimilarity matrix computed from pairwise alignments
- a structure tree inferred from the SDM matrix
- sequence extraction directly from the input structures
- `MAFFT` for multiple sequence alignment
- `IQ-TREE 3` for maximum-likelihood sequence phylogeny
- summary QC, distance-matrix comparison, and heatmap/report outputs

## Environment

The repository includes a `mamba` environment specification in [`environment.yml`](./environment.yml).

Create and activate it with:

```bash
mamba env create -f environment.yml
mamba activate structphylogeny
```

## Input

Place protein structures in a directory such as `data/my_structures/`.

Supported structure formats:

- `.cif`
- `.mmcif`
- `.pdb`
- `.ent`

The current default dataset is:

```text
data/mouse_lipocalin_structures
```

## Run

The default configuration is in [`config/config.yaml`](./config/config.yaml).

Run the full pipeline with:

```bash
snakemake --cores 8
```

Or dry-run it first:

```bash
snakemake --cores 8 --dry-run
```

## Outputs

Main outputs land under `results/`:

- `results/manifest/structures.tsv`: input manifest and basic QC
- `results/sequences/proteins.fasta`: extracted sequences
- `results/sequences/proteins.aligned.fasta`: MAFFT alignment
- `results/trees/sequence.treefile`: IQ-TREE sequence phylogeny
- `results/structure/gtalign_pairwise.tsv`: parsed GTalign pairwise metrics
- `results/structure/sdm.tsv`: structural dissimilarity matrix
- `results/trees/structure_sdm.nwk`: structure tree from the SDM matrix
- `results/reports/sdm_heatmap.png`: SDM heatmap
- `results/reports/sequence_distance_heatmap.png`: sequence p-distance heatmap
- `results/reports/matrix_correlation.tsv`: agreement between structure and sequence distances
- `results/reports/report.md`: concise summary report

## Notes

- The paper used DALI and a PHYLIP distance-tree method. This repository uses `GTalign` for structural comparison, as it more modern, faster and better tool.
- The SDM implementation follows the Lakshmi et al. formula, but approximates PFTE from GTalign's aligned-residue count and the smaller structure length.
- `IQ-TREE 3` is used for the sequence phylogeny. The structure phylogeny is inferred from the SDM distance matrix with a neighbor-joining distance-tree method.
- By default, if a structure has multiple chains, the pipeline uses the longest protein chain.

## Tests

Run the lightweight unit tests with:

```bash
PYTHONPATH=src pytest
```
