# protein-design-ops

Backbone-conditioned sequence design with **independent rescoring**: given a
protein backbone, ProteinMPNN generates candidate sequences, then ESM-2
scores each candidate's zero-shot sequence fitness, a second model's
opinion on the first model's output. Candidates both models agree on are
the defensible shortlist. Disagreement is reported, not hidden.

## Pipeline

```
backbone -> generate (ProteinMPNN) -> score (ESM-2) -> report
```

- `backbone.py` - parse PDB, extract design chain + native sequence
- `generate.py` - subprocess driver for upstream ProteinMPNN (external
  clone, not vendored, see `config.mpnn.repo_path`), parses FASTA headers
  carrying MPNN's own score and sequence recovery
- `score.py` - ESM-2 mean pseudo-log-likelihood per candidate (mask each
  position, log-prob of the actual residue, average)
- `report.py` - consensus ranking, score correlation, diversity,
  identity-to-WT, scatter figure

## Quickstart

```bash
git clone https://github.com/dauparas/ProteinMPNN ../proteinmpnn-ext
uv venv --python 3.11 && uv pip install -e ".[dev]"
.venv/bin/python -m snakemake --cores 2
```

Demo backbone: trp-cage miniprotein (1L2Y, 20 aa). The whole DAG runs in
seconds on CPU, including ESM-2 scoring. A second config
(`config/config_1ubq.yaml`, ubiquitin, 76 aa natural fold) runs with
`DESIGN_CONFIG=config/config_1ubq.yaml snakemake`. Outputs are keyed by
`backbone.id` under `results/<id>/`.

## Observed results (16 candidates @ T=0.1 each)

| Backbone | Fold | Spearman(MPNN, ESM-2) | Mean recovery | Native ESM-2 | Designed ESM-2 mean |
|---|---|---:|---:|---:|---:|
| 1L2Y trp-cage (20 aa) | engineered miniprotein | **-0.632** | 0.38 | -2.87 | -2.59 |
| 1UBQ ubiquitin (76 aa) | natural globular | **-0.529** | 0.55 | -2.32 | -1.92 |

- On both backbones the two models **anticorrelate**, so the consensus shortlist is a
  compromise, not the best on either metric. This is reported as-is.
  Challenging single-model confidence is why the second scorer exists.
- Designed sequences out-score the native on ESM-2 for both backbones.
  For the natural fold the gap is wider (-1.92 vs -2.32). Plausible for
  fixed-backbone redesign. It is a sequence-fitness observation, not a
  folding or function claim.

`design_report.json` carries a `provenance` block: backbone id,
ProteinMPNN upstream commit + weights, sampling params/seed, and the
ESM-2 model id, enough to reproduce a run exactly.

**Reproducibility caveat:** ProteinMPNN's `--seed`
treats `0` as "pick a random seed" (the upstream script does
`if args.seed:`, and 0 is falsy). An earlier config set `seed: 0` believing
it pinned sampling. The runs above use `seed: 37` and are deterministic
(bit-identical candidates across repeated runs). The earlier ubiquitin
correlation near zero was a random-seed artifact, not a real
backbone-dependence signal. With sampling pinned, both backbones
anticorrelate.

## Scope honesty

- Candidates are computational designs only. No claim of stability,
  folding, or function. That requires structure prediction of the
  designs (e.g. ESMFold/Boltz) or experiment.
- Two backbones, one temperature, one seed is a demo, not a study.
  `config/config.yaml` holds all settings for sweeps.
- ProteinMPNN and ESM-2 are upstream models (dauparas/ProteinMPNN,
  facebook/esm). This repo is orchestration + evaluation.

MIT licensed.
