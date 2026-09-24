"""Sequence generation: invoke upstream ProteinMPNN and parse its output.

ProteinMPNN is external (config.mpnn.repo_path); we call its
protein_mpnn_run.py in a subprocess and parse the emitted FASTA, which
carries the model's own score and seq_recovery in each header.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

CONFIG = os.environ.get("DESIGN_CONFIG", "config/config.yaml")


def parse_mpnn_fasta(path: str) -> list[dict]:
    """ProteinMPNN FASTA -> records; first entry is the native sequence."""
    records = []
    header = None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
        elif line and header is not None:
            meta = dict(re.findall(r"(\w+)=(\S+?)(?:,|$)", header))
            records.append(
                {
                    "header": header,
                    "seq": line,
                    "mpnn_score": float(meta["score"]),
                    "seq_recovery": float(meta.get("seq_recovery", 1.0)),
                    "is_native": "designed_chains" in header
                        and "sample" not in meta,
                }
            )
            header = None
    return records


def generate(cfg: dict, backbone: dict, workdir: Path) -> list[dict]:
    m = cfg["mpnn"]
    out_dir = workdir / "mpnn_out"
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_abs = str(Path(backbone["pdb"]).resolve())
    cmd = [
        sys.executable,
        str(Path(m["repo_path"]) / m["run_script"]),
        "--pdb_path", pdb_abs,
        "--pdb_path_chains", backbone["chain"],
        "--out_folder", str(out_dir),
        "--num_seq_per_target", str(m["num_seq_per_target"]),
        "--sampling_temp", str(m["sampling_temp"]),
        "--batch_size", str(m["batch_size"]),
        "--seed", str(m["seed"]),
        "--model_name", m["model_name"],
        "--path_to_model_weights",
        str(Path(m["repo_path"]) / Path(m["weights"]).parent) + "/",
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    fasta = next((out_dir / "seqs").glob("*.fa"))
    return parse_mpnn_fasta(str(fasta))


def main() -> None:
    backbone_path, out_path = sys.argv[1], sys.argv[2]
    with open(CONFIG) as f:
        cfg = yaml.safe_load(f)
    backbone = json.load(open(backbone_path))
    with tempfile.TemporaryDirectory() as tmp:
        records = generate(cfg, backbone, Path(tmp))
    with open(out_path, "w") as f:
        json.dump(records, f, indent=2)
    n_native = sum(r["is_native"] for r in records)
    print(
        f"generate: {len(records)} sequences "
        f"({len(records) - n_native} designed + {n_native} native) "
        f"-> {out_path}"
    )


if __name__ == "__main__":
    main()
