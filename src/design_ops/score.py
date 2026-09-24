"""Independent rescoring: ESM-2 mean pseudo-log-likelihood per sequence.

For each candidate, mask every position in turn and take the model's
log-probability of the actual residue; the mean over positions is the
sequence's zero-shot fitness score. This is an independent second opinion
on MPNN's own score — the cross-model agreement is the report's signal.
"""

from __future__ import annotations

import json
import os
import sys

import yaml

CONFIG = os.environ.get("DESIGN_CONFIG", "config/config.yaml")


def pll(seq: str, tok, model) -> float:
    """Mean masked-marginal log-prob of a sequence. Pure w.r.t. model+tok."""
    import torch

    total = 0.0
    for i in range(len(seq)):
        masked = seq[:i] + tok.mask_token + seq[i + 1 :]
        inputs = tok(masked, return_tensors="pt")
        with torch.no_grad():
            logits = model(**inputs).logits[0]
        mask_idx = (
            inputs["input_ids"][0] == tok.mask_token_id
        ).nonzero()[0].item()
        lp = torch.log_softmax(logits[mask_idx], dim=-1)
        total += float(lp[tok.convert_tokens_to_ids(seq[i])])
    return total / len(seq)


def score_candidates(records: list[dict], model_name: str) -> list[dict]:
    from transformers import AutoTokenizer, EsmForMaskedLM

    tok = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name).eval()
    for r in records:
        r["esm_pll"] = round(pll(r["seq"], tok, model), 4)
    return records


def main() -> None:
    in_path, out_path = sys.argv[1], sys.argv[2]
    with open(CONFIG) as f:
        cfg = yaml.safe_load(f)
    records = json.load(open(in_path))
    records = score_candidates(records, cfg["esm"]["model"])
    with open(out_path, "w") as f:
        json.dump(records, f, indent=2)
    designed = [r for r in records if not r["is_native"]]
    print(
        f"score: {len(designed)} designed + "
        f"{len(records) - len(designed)} native scored -> {out_path}"
    )


if __name__ == "__main__":
    main()
