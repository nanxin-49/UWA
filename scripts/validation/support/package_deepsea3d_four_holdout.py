"""Package only the frozen inputs needed for the four-frequency remote run."""
import hashlib
import json
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_four_holdout_20261026"
OLD = ROOT / "results/validation/deepsea3d_5k7k_20261026/input_manifest.json"
NEW = [
    "scripts/validation/validate_deepsea3d_four_holdout_frequency.m",
    "scripts/validation/run_deepsea3d_four_holdout_frequency.m",
    "results/validation/deepsea3d_four_holdout_20261026/prediction_lock.json",
    "results/validation/deepsea3d_four_holdout_20261026/predictions_pre_FE.csv",
    "results/validation/deepsea3d_four_holdout_20261026/launch_ubuntu_matlab.sh",
]


def main():
    rows = json.loads(OLD.read_text(encoding="utf-8-sig"))
    rows = [row for row in rows if row["path"] not in {
        "scripts/validation/validate_deepsea3d_5k7k_frequency.m",
        "scripts/validation/run_deepsea3d_5k7k_frequency.m",
        "scripts/validation/check_deepsea3d_ubuntu_compatibility.m",
    }]
    for relative in NEW:
        path = ROOT / relative
        rows.append({"path": relative, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                     "bytes": path.stat().st_size})
    for row in rows:
        path = ROOT / row["path"]
        assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], row["path"]
    archive = OUT / "remote_inputs.tar"
    assert not archive.exists(), "Never replace a transferred input archive"
    with tarfile.open(archive, "w") as tar:
        for row in rows:
            tar.add(ROOT / row["path"], arcname=row["path"], recursive=False)
    (OUT / "remote_input_manifest.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    (OUT / "remote_input.sha256").write_bytes("".join(
        f'{row["sha256"]}  {row["path"]}\n' for row in rows).encode("utf-8"))
    print(json.dumps({"input_count": len(rows), "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
                      "archive_bytes": archive.stat().st_size}))


if __name__ == "__main__":
    main()
