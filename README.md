# protein-design-ops

Backbone-conditioned sequence design with **independent rescoring**: given a
protein backbone, ProteinMPNN generates candidate sequences, then ESM-2
scores each candidate's zero-shot sequence fitness — a second model's
opinion on the first model's output. Candidates both models agree on are
the defensible shortlist; disagreement is reported, not hidden.

## Pipeline

```
backbone -> generate (ProteinMPNN) -> score (ESM-2) -> report
```

- `backbone.py` — parse PDB, extract design chain + native sequence
- `generate.py` — subprocess driver for upstream ProteinMPNN (external
  clone, not vendored — see `config.mpnn.repo_path`), parses FASTA headers
  carrying MPNN's own score and sequence recovery
- `score.py` — ESM-2 mean pseudo-log-likelihood per candidate: mask each
  position, log-prob of the actual residue, average
- `report.py` — consensus ranking, score correlation, diversity,
  identity-to-WT, scatter figure

## Quickstart

```bash
git clone https://github.com/dauparas/ProteinMPNN ../proteinmpnn-ext
uv venv --python 3.11 && uv pip install -e ".[dev]"
.venv/bin/python -m snakemake --cores 2
```

Demo backbone: trp-cage miniprotein (1L2Y, 20 aa) — the whole DAG runs in
seconds on CPU, including ESM-2 scoring.

## Observed result (1L2Y, 16 candidates @ T=0.1)

- Mean pairwise identity 0.71, mean sequence recovery 37.5%
- **MPNN and ESM-2 scores anticorrelate (spearman -0.52)** — on this
  backbone the models genuinely disagree; the consensus shortlist is a
  compromise, not the best on either metric. This is reported rather than
  smoothed over — single-model confidence is exactly what the second
  scorer exists to challenge.
- Native trp-cage scores *worse* on ESM-2 (-2.87) than the designed mean
  (-2.65) — consistent with it being an engineered miniprotein rather than
  a natural fold.

## Scope honesty

- Candidates are computational designs only. No claim of stability,
  folding, or function — that requires structure prediction of the
  designs (e.g. ESMFold/Boltz) or experiment.
- One backbone, one temperature, one seed is a demo, not a study.
  `config/config.yaml` holds all settings for sweeps.
- ProteinMPNN and ESM-2 are upstream models (dauparas/ProteinMPNN,
  facebook/esm); this repo is orchestration + evaluation.

MIT licensed.
