"""Verify and extract only the task-specific returned result bundle."""
import hashlib
import io
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_65grid_four_20261027"
ARCHIVE = OUT / "transfer_bundle.tar"
EXPECTED = "09ac7f0ece0b7aa7617a2e0350d7ca7037f370857d848780cf132ea5593ac653"
FIGURE_EXPECTED = "d676ad14f8e6acb71089ca152e8d45795b0da00f315afc8acb093d8abb25c537"
FIGURE_V3_EXPECTED = "f2bc9993d42c4d81fd92a02acc80a48b04a9520441b3ef7c4a4c8c83b8b35dec"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def ingest(archive):
    count = 0
    with tarfile.open(archive, "r") as handle:
        members = handle.getmembers()
        for member in members:
            name = member.name.removeprefix("./")
            target = OUT / name
            assert name and ".." not in Path(name).parts and not Path(name).is_absolute()
            assert target.resolve().is_relative_to(OUT.resolve())
            assert member.isdir() or member.isfile(), name
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                data = handle.extractfile(member).read()
                if target.exists():
                    assert target.read_bytes() == data, name
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with target.open("xb") as stream:
                        stream.write(data)
                count += 1
    return count

def main():
    assert sha(ARCHIVE) == EXPECTED
    count = ingest(ARCHIVE)
    for line in (OUT / "output.sha256").read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        assert sha(OUT / relative) == digest, relative
    figure_bundle = OUT / "figure_bundle.tar"
    figure_count = 0
    if figure_bundle.exists():
        assert sha(figure_bundle) == FIGURE_EXPECTED
        figure_count = ingest(figure_bundle)
        for line in (OUT / "figure_output.sha256").read_text(encoding="utf-8").splitlines():
            digest, relative = line.split("  ", 1)
            assert sha(OUT / relative) == digest, relative
    final_bundle = OUT / "figure_bundle_v3.tar"
    final_count = 0
    if final_bundle.exists():
        assert sha(final_bundle) == FIGURE_V3_EXPECTED
        final_count = ingest(final_bundle)
        for line in (OUT / "figure_output_v3.sha256").read_text(encoding="utf-8").splitlines():
            digest, relative = line.split("  ", 1)
            assert sha(OUT / relative) == digest, relative
    print({"bundle_sha256": EXPECTED, "returned_files_verified": len((OUT / "output.sha256").read_text(encoding="utf-8").splitlines()),
           "extracted_file_count": count, "figure_bundle_verified": bool(figure_count),
           "figure_files_extracted": figure_count, "final_figure_bundle_verified": bool(final_count)})

if __name__ == "__main__":
    main()
