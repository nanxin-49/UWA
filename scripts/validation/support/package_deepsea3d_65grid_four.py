"""Package only SHA-verified inputs for the isolated Ubuntu four-point run."""
import hashlib
import json
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_65grid_four_20261027"
PRIOR = ROOT / "results/validation/deepsea3d_four_holdout_20261026/remote_input_manifest.json"
ADDED = [
    "scripts/validation/validate_deepsea3d_65grid_four_frequency.m",
    "scripts/validation/run_deepsea3d_65grid_four_frequency.m",
    "scripts/validation/support/launch_deepsea3d_65grid_four_remote.sh",
    "results/validation/deepsea3d_65grid_four_20261027/prediction_lock.json",
    "results/validation/deepsea3d_65grid_four_20261027/predictions_pre_FE.csv",
    "results/validation/deepsea3d_65grid_four_20261027/prediction_lock_audit.json",
]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    rows = json.loads(PRIOR.read_text(encoding="utf-8-sig"))[:14]
    for row in rows:
        p = ROOT / row["path"]
        assert p.is_file() and p.stat().st_size == row["bytes"] and sha(p) == row["sha256"], row["path"]
    for name in ADDED:
        p = ROOT / name
        assert p.is_file()
        rows.append({"path": name, "sha256": sha(p), "bytes": p.stat().st_size})
    archive = OUT / "remote_inputs.tar"
    assert not archive.exists()
    with tarfile.open(archive, "w") as handle:
        for row in rows:
            handle.add(ROOT / row["path"], arcname=row["path"], recursive=False)
    (OUT / "remote_input_manifest.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    (OUT / "remote_input.sha256").write_bytes("".join(
        f'{row["sha256"]}  {row["path"]}\n' for row in rows).encode("ascii"))
    print(json.dumps({"file_count": len(rows), "archive_sha256": sha(archive),
                      "archive_bytes": archive.stat().st_size}))

if __name__ == "__main__":
    main()
